# This file takes as input a filename (in the variable fname) that
# describes a circuit in netlist format (similar to SPICE), and then
# performs a symbolic analysis of the circuit using Modified Nodal Analysis
# (MNA).  A full description of MNA and how to use this code is at:
# http://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA1.html
#
# Written by Erik Cheever, Swarthmore College, 2019
# echeeve1@swarthmore.edu
#
# https://github.com/echeever/scam
#
# Julia version by Jakub Ladman ladmanj@volny.cz, 2022

using Printf
using TickTock
using Symbolics

function parse_string_matrix(M)

    M = Meta.parse.(M)
    M = parse_expr_to_symbolic.(M, (Main,))
    M = convert(Matrix{Num}, M)
    M = simplify.(M)

end

"""
    solve_circuit(fname)

TBW
"""
function solve_circuit(fname)

    @printf("\nStarted -- please be patient.\n\n")

    #fname = "ex.cir"
    #fname = "circuits/thermal2_sym.cir"
    # Print out the netlist (a file describing the circuit with one circuit
    # per line.
    @printf("Netlist: %s\n", fname)

    fid = open(fname)
    fileIn = split.(readlines(fid))  # Read file (up to 6 items per line
    # Split each line into 6 columns, the meaning of the last 3 columns will
    # vary from element to element.  The first item is always the name of the
    # circuit element and the second and third items are always node numbers.

    display(fileIn)

    Name = map(x -> get(x, 1, ""), fileIn)
    N1 = map(x -> get(x, 2, ""), fileIn)
    N2 = map(x -> get(x, 3, ""), fileIn)
    arg3 = map(x -> get(x, 4, ""), fileIn)
    arg4 = map(x -> get(x, 5, ""), fileIn)
    arg5 = map(x -> get(x, 6, ""), fileIn)
    # Name, node1, node2, and up to three other arguments.
    close(fid)

    nLines = length(Name)  # Number of lines in file (or elements in circuit).

    N1 = parse.(Int, N1)   # Get node numbers
    N2 = parse.(Int, N2)

    tick()                  # Begin timing.

    n = maximum(abs, [N1; N2], dims=1)[1]   # Find highest node number (i.e., number of nodes)

    m = 0 # "m" is the number of voltage sources, determined below.

    for k1 in 1:nLines                  # Check all lines to find voltage sources
        if ((Name[k1][1] == 'V')
            || (Name[k1][1] == 'O')
            || (Name[k1][1] == 'E')
            || (Name[k1][1] == 'H'))
            # These are the circuit elements with
            m = m + 1             # We have voltage source, increment m.
        end
    end

    # Preallocate all arrays (use Litovski's notation).
    G = string.(zeros(Int, n, n))    # G is nxn filled with '0'
    B = string.(zeros(Int, n, m))
    C = string.(zeros(Int, m, n))
    D = string.(zeros(Int, m, m))
    i = string.(zeros(Int, n, 1))
    e = string.(zeros(Int, m, 1))
    j = string.(zeros(Int, m, 1))
    v = string.("v_", 1:n)          # v is filled with node names

    # We need to keep track of the number of voltage sources we've parsed
    # so far as we go through file.  We start with zero.
    vsCnt = 0

    # This loop does the bulk of filling in the arrays.  It scans line by line
    # and fills in the arrays depending on the type of element found on the
    # current line.
    # See http://lpsa.swarthmore.edu/Systems/Electrical/mna/MNA3.html for details.
    for k1 in 1:nLines
        n1 = N1[k1]   # Get the two node numbers
        n2 = N2[k1]

        if ((Name[k1][1] == 'R')
            || (Name[k1][1] == 'L')
            || (Name[k1][1] == 'C'))
            # Passive element
            # RXXXXXXX N1 N2 VALUE
            if (Name[k1][1] == 'R')  # Find 1/impedance for each element type.
                g = string.("1/", Name[k1])
            elseif (Name[k1][1] == 'L')
                g = string.("1/s/", Name[k1])
            elseif (Name[k1][1] == 'C')
                g = string.("s*", Name[k1])
            end
            # Here we fill in G array by adding conductance.
            # The procedure is slightly different if one of the nodes is
            # ground, so check for those accordingly.
            if (n1 == 0)
                G[n2, n2] = string.(G[n2, n2], " + ", g)  # Add conductance.
            elseif (n2 == 0)
                G[n1, n1] = string.(G[n1, n1], " + ", g)  # Add conductance.
            else
                G[n1, n1] = string.(G[n1, n1], " + ", g)  # Add conductance.
                G[n2, n2] = string.(G[n2, n2], " + ", g)  # Add conductance.
                G[n1, n2] = string.(G[n1, n2], " - ", g)  # Sub conductance.
                G[n2, n1] = string.(G[n2, n1], " - ", g)  # Sub conductance.
            end

            # Independent voltage source.
        elseif (Name[k1][1] == 'V') # VXXXXXXX N1 N2 VALUE    (N1=anode, N2=cathode)
             vsCnt = vsCnt + 1  # Keep track of which source this is.
            # Now fill in B and C arrays (again, process is slightly
            # different if one of the nodes is ground).
            if (n1 != 0)
                B[n1, vsCnt] = string.(B[n1, vsCnt], " + 1")
                C[vsCnt, n1] = string.(C[vsCnt, n1], " + 1")
            end
            if (n2 != 0)
                B[n2, vsCnt] = string.(B[n2, vsCnt], " - 1")
                C[vsCnt, n2] = string.(C[vsCnt, n2], " - 1")
            end
            e[vsCnt] = Name[k1]         # Add Name of source to RHS
            j[vsCnt] = string.("I_", Name[k1])  # Add current through source to unknowns


        # Independent current source
        elseif (Name[k1][1] == 'I') # IXXXXXXX N1 N2 VALUE  (Current N1 to N2)
            # Add current to nodes (if node is not ground)
            if (n1 != 0)
                i[n1] = string.(i[n1], " - ", Name[k1]) # subtract current from n1
            end
            if (n2 != 0)
                i[n2] = string.(i[n2], " + ", Name[k1]) # add current to n2
            end

            # Op amp
        elseif (Name[k1][1] == 'O')  # 0XXXXXXX N1 N2 N3 VALUE  (N1=+, N2=-, N3=Vout)
            n3 = parse.(Int, arg3[k1])  # This find n3
             vsCnt = vsCnt + 1          # Keep track of number of voltage sources

            # Change B and C matrices as appropriate.
            B[n3, vsCnt] = string.(B[n3, vsCnt], " + 1")
            if (n1 != 0)
                C[vsCnt, n1] = string.(C[vsCnt, n1], " + 1")
            end
            if (n2 != 0)
                C[vsCnt, n2] = string.(C[vsCnt, n2], " - 1")
            end
            j[vsCnt] = string.("I_", Name[k1])  # Add current through source to unknowns


        # Voltage controlled voltage source
        elseif (Name[k1][1] == 'E') # VCVS
             vsCnt = vsCnt + 1           # Keep track of number of voltage sources
            nc1 = parse.(Int, arg3[k1])  # Control voltage, pos side
            nc2 = parse.(Int, arg4[k1])  # Control voltage, neg side

            # Change B and C matrices as appropriate for output nodes.
            #  (if node is not ground)
            if (n1 != 0)
                B[n1, vsCnt] = string.(B[n1, vsCnt], " + 1")
                C[vsCnt, n1] = string.(C[vsCnt, n1], " + 1")
            end
            if (n2 != 0)
                B[n2, vsCnt] = string.(B[n2, vsCnt], " - 1")
                C[vsCnt, n2] = string.(C[vsCnt, n2], " - 1")
            end

            # Change C matrix as appropriate for input nodes
            # (if node is not ground)
            if (nc1 != 0)
                C[vsCnt, nc1] = string.(C[vsCnt, nc1], " - ", Name[k1])
            end
            if (nc2 != 0)
                C[vsCnt, nc2] = string.(C[vsCnt, nc2], " + ", Name[k1])
            end

            j[vsCnt] = string.("I_", Name[k1]) # Add current through source to unknowns


        # Voltage controlled current source
        elseif (Name[k1][1] == 'G')    # VCCS GXXXXXXX N+ N- NC+ NC- VALUE
            nc1 = parse.(Int, arg3[k1])    # Control voltage, pos side
            nc2 = parse.(Int, arg4[k1])    # Control voltage, neg side
            g = Name[k1]

            # Create a string that shows if each node is ~= zero
            # (i.e., we find which nodes are grounded).
            myString = string.((
                    if (n1 != 0)
                        "1"
                    else
                        "0"
                    end
                ), (
                    if (n2 != 0)
                        "1"
                    else
                        "0"
                    end
                ), (
                    if (nc1 != 0)
                        "1"
                    else
                        "0"
                    end
                ), (
                    if (nc2 != 0)
                        "1"
                    else
                        "0"
                    end
                ))
            #  myString = myString(~isspace(myString))  # Remove spaces
            # Checking all different conditions gets complicated.  There
            # may be a simpler way, but here we just brute force it and
            # check all 16 possible.
            if ((myString == "0000")
                || (myString == "0011")
                || (myString == "0001")
                || (myString == "0010")
                || (myString == "0100")
                || (myString == "1000")
                || (myString == "1100"))
                error("error in VCCS") # This should never happen
            end

            if (myString == "1111")  # All nodes are non-zero
                G[n1, nc1] = string.(G[n1, nc1], " + ", g)
                G[n1, nc2] = string.(G[n1, nc2], " - ", g)
                G[n2, nc1] = string.(G[n2, nc1], " - ", g)
                G[n2, nc2] = string.(G[n2, nc2], " + ", g)
            elseif (myString == "0111")  # n1 is zero - so don't include
                G[n2, nc1] = string.(G[n2, nc1], " - ", g)
                G[n2, nc2] = string.(G[n2, nc2], " + ", g)
            elseif (myString == "0101")
                G[n2, nc2] = string.(G[n2, nc2], " + ", g)
            elseif (myString == "0110")
                G[n2, nc1] = string.(G[n2, nc1], " - ", g)
            elseif (myString == "1011")
                G[n1, nc1] = string.(G[n1, nc1], " + ", g)
                G[n1, nc2] = string.(G[n1, nc2], " - ", g)
            elseif (myString == "1001")
                G[n1, nc2] = string.(G[n1, nc2], " - ", g)
            elseif (myString == "1010")
                G[n1, nc1] = string.(G[n1, nc1], " + ", g)
            elseif (myString == "1101")
                G[n1, nc2] = string.(G[n1, nc2], " - ", g)
                G[n2, nc2] = string.(G[n2, nc2], " + ", g)
            elseif (myString == "1110")
                G[n1, nc1] = string.(G[n2, nc1], " + ", g)
                G[n2, nc1] = string.(G[n2, nc1], " - ", g)
            end

            # Current controlled current source.
            # For the CCCS we need the controlling current, which is
            # defined as the current through one of the voltage sources.
            # Since this voltage may not have been defined yet (i.e., it
            # comes later in the circuit definition file), we leave this
            # part of the matrix generation for later.
            # For the CCCS there is nothing to add at this point.
        elseif (Name[k1][1] == 'F')    # CCCS FXXXXXXX N+ N- VNAM VALUE
        # Current controlled voltage source
        # For the CCVS we need the controlling current which is defined as the
        # current through one of the voltage sources.  Since this voltage may not
        # have been defined yet (i.e., it comes later in the circuit definition
        # file), we leave this part of the matrix generation for later.
        elseif (Name[k1][1] == 'H')    # CCVS
             vsCnt = vsCnt + 1 # Keep track of number of voltage sources
            # Change B and C as appropriate (if node is not ground)
            if (n1 != 0)
                B[n1, vsCnt] = string.(B[n1, vsCnt], " + 1")
                C[vsCnt, n1] = string.(C[vsCnt, n1], " + 1")
            end
            if (n2 != 0)
                B[n2, vsCnt] = string.(B[n2, vsCnt], " - 1")
                C[vsCnt, n2] = string.(C[vsCnt, n2], " - 1")
            end
            j[vsCnt] = string.("I_", Name[k1]) # Add current through source to unknowns

        end
    end

    # At this point all voltage sources have been parsed.  We can now go
    # through and finish off the CCVS and CCCS elements (which depend on the
    # current through those sources).
    for k1 in 1:nLines
        n1 = N1[k1]
        n2 = N2[k1]
        if (Name[k1][1] == 'H')
            # Here we find the indices in the matrix j:
            #    of the controlling voltage (cvInd)
            #    as well as the index of this element (hInd)
            cv = arg3[k1]  # Name of controlling voltages
            cvInd = findall(contains(j, cv))  # Index of controlling voltage.
            hInd = findall(contains(j, Name[k1])) # Index of CCVS (this element)
            D[hInd, cvInd] = string.('-', Name[k1])  # Set the value of the D matrix.
        elseif (Name[k1][1] == 'F')
            # Here we find the index in the matrix j of the controlling
            # voltage (cvInd)
            cv = arg3[k1] # Name of controlling voltages
            cvInd = findall(contains(j, cv))  # Index of controlling voltage
            # Set the B matrix accordingly.
            if (n1 != 0)
                B[n1, cvInd] = string.(B[n1, cvInd], " + ", Name[k1])
            end
            if (n2 != 0)
                B[n2, cvInd] = string.(B[n2, cvInd], " - ", Name[k1])
            end
        end
    end
    ##  The submatrices are now complete.  Form the A, x, and z matrices,
    # and solve!

    A = parse_string_matrix([G B; C D]) # Create and display A matrix
    @printf("\nThe A matrix: \n")
    display(A)

    x = parse_string_matrix([v; j]) # Create and display x matrix
    @printf("\nThe x matrix: \n")
    display(x)

    z = parse_string_matrix([i; e]) # Create and display z matrix
    @printf("\nThe z matrix:  \n")
    display(z)

    # Find all variables in matrices (symvar) and make them symbolic (syms)
    #syms([symvar(A), symvar(x), symvar(z)])

    # Display the matrix equation
    @printf("\nThe matrix equation: \n")
    display(A * x ~ z)

    a = Symbolics.simplify.(A \ z)  # Get the solution, this is the heart of the algorithm.

    @show a
    #render(latexify(a))

    # instead of polluting variable scope, let's put eqn system solution into a hash
    #for i in 1:length(a)  # Assign each solution to its output variable.
    #    eval(sprintf("%s = %s;",x(i),a(i)))
    #end
    y = Dict(Symbol.(x) .=> a)
    # instead of working with e.g. v_1 we use y[:v_1]

    @printf("\nThe solution:  \n")
    @show y

    # Lastly, assign any numeric values to symbolic variables.
    # Go through the elements a line at a time and see if the values are valid
    # numbers.  If so, assign them to the variable name.  Then you can use
    # "eval" to get numerical results.
    for k1 in 1:nLines
        num = nothing
        status = false
        if ((Name[k1][1] == 'R')
            || (Name[k1][1] == 'L')
            || (Name[k1][1] == 'C')
            || (Name[k1][1] == 'V')
            || (Name[k1][1] == 'I'))
            # These circuit elements defined by three variables, 2 nodes and a
            # value.  The third variable (arg3) is the value.

            status = try
                num = parse.(Float64, arg3[k1]) # ok<ST2NM>
                status = true
            catch
                status = false
            end
            # Elements defined by four variables, arg4 is the value.
        elseif ((Name[k1][1] == 'H')
                ||
                (Name[k1][1] == 'F'))

            status = try
                num = parse.(Float64, arg4[k1]) # ok<ST2NM>
                status = true
            catch
                status = false
            end
            # Elements defined by five variables, arg5 is the value.
        elseif ((Name[k1][1] == 'E')
                ||
                (Name[k1][1] == 'G'))

            status = try
                num = parse.(Float64, arg3[k1]) # ok<ST2NM>
                status = true
            catch
                status = false
            end
        end
        #@show status num
        if status  # status will be true if the argument was a valid number.
            # If the number is valid, assign it to the variable.
            #eval(string.(Name[k1], " = ", num))
            # let's substitute into the eqn solution expressions instead:
            var = parse_expr_to_symbolic.(Symbol(Name[k1]), (Main,))[1]
            for ex in keys(y)
                y[ex] = SymbolicUtils.substitute(y[ex], Dict([var => num]))
            end
        end
    end

    @printf("\nElapsed time is %g seconds.\n", tok())
end

# solve_circuit("circuits/thermal2_sym.cir")
# solve_circuit("circuits/example0.cir")
# solve_circuit("circuits/example1.cir")
# solve_circuit("circuits/example3.cir")
# solve_circuit("circuits/example4.cir")
solve_circuit("circuits/example5.cir")
# solve_circuit("circuits/example6.cir")
# solve_circuit("circuits/example7.cir")
# solve_circuit("circuits/exampleE.cir")
# solve_circuit("circuits/exampleF.cir")
# solve_circuit("circuits/exampleG.cir")
# solve_circuit("circuits/exampleH.cir")

