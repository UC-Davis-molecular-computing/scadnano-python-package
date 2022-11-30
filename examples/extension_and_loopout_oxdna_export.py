import scadnano as sc


def main() -> None:
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
    d.write_oxdna_files(directory='oxdna')


def create_design() -> sc.Design:
    width = 21
    helices = [sc.Helix(max_offset=32) for _ in range(3)]
    helices[0].roll = 210
    helices[1].roll = 20
    helices[2].roll = 210
    design = sc.Design(helices=helices, grid=sc.square)

    design.draw_strand(0, 0).extension_5p(5) \
        .move(width) \
        .cross(1).move(-width) \
        .loopout(2, 32) \
        .move(width) \
        .extension_3p(7)

    return design


if __name__ == '__main__':
    main()
