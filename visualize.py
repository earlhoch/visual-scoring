from pymol import cmd, stored

#Colors the provided ligand and receptor molecules according to their
#b-factor values
def visualize( lig, rec ):


    stored.b_factors = []

    cmd.iterate_state(1, selector.process(lig), "stored.b_factors.append(b)")
    cmd.iterate_state(1, selector.process(rec), "stored.b_factors.append(b)")

    max = 0.00

    for value in stored.b_factors:
        if abs(value) > max:
            max = abs(value)

    print("Maximum absolute b-factor value: %s") % (max)

    #Uses same magnitude for maximum and minimum to stay symmetrical around zero
    cmd.spectrum("b", "red_white_green", minimum = (0-max), maximum = max, selection = lig)
    cmd.spectrum("b", "red_white_green", minimum = (0-max), maximum = max, selection = rec)


cmd.extend( "visualize", visualize );
