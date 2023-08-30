import scadnano as sc


def create_design() -> sc.Design:
    # shows how to make consecutive domains on a helix, separated by a crossover that appears horizontal
    # this is useful when doing single-stranded tile designs, for instance, or any other design
    # where we have consecutive domains on a single helix.
    #
    # 0 [------+^+------>
    #
    # 1 <------+^+------]
    design = sc.Design(helices=[sc.Helix(100), sc.Helix(100)], grid=sc.square)
    design.draw_strand(0, 0).move(8).move(8)
    design.draw_strand(1, 16).move(-8).move(-8)

    # XXX: the following code raises an exception because it tries to add a crossover where
    # there already is one. This can be surprising since it would work with the following
    # similar-looking design that has a single longer domain per strand
    #
    # 0 [------------>
    #
    # 1 <------------]
    # design.add_full_crossover(helix=0, helix2=1, offset=8, forward=True)

    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
