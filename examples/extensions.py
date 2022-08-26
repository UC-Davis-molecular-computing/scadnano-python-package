import scadnano as sc


def main() -> None:
    d = create_design()
    d.write_scadnano_file(directory='output_designs')


def create_design() -> sc.Design:
    width = 8
    helices = [sc.Helix(max_offset=32) for _ in range(3)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.draw_strand(0, 0).extension_5p(5, display_length=2.5, display_angle=45) \
        .move(width).cross(1).move(-width).loopout(2, 3).move(width) \
        .extension_3p(7).with_domain_name("ext_3p")

    design.draw_strand(0, 24).extension_5p(5, display_length=3.5, display_angle=60) \
        .move(-width).cross(1).move(width).loopout(2, 3).move(-width) \
        .extension_3p(7).with_domain_name("ext_3p_top")

    return design


if __name__ == '__main__':
    main()
