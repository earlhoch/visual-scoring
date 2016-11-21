from pymol import cmd, stored

#Colors the provided ligand and receptor molecules according to their
#b-factor values
def visualize( *args, **kwargs ):

    if len(args) < 1:
           print("No molecules selected")
           return



    stored.b_factors = []
    for selection in args:

            cmd.iterate_state(1, selector.process(selection), "stored.b_factors.append(b)")

    max = 0.00

    for value in stored.b_factors:
        if abs(value) > max:
            max = abs(value)

    print("Maximum absolute b-factor value: %s") % (max)

    #Uses same magnitude for maximum and minimum to stay symmetrical around zero
    for selection_item in args:

            cmd.spectrum("b", "red_white_green", minimum = (0-max), maximum =
                            max, selection = selection_item)


cmd.extend( "visualize", visualize );
