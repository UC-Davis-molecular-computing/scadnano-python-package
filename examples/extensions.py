import scadnano as sc


def main() -> None:
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
    d.write_oxdna_files(directory='oxdna')


def create_design() -> sc.Design:
    width = 8
    helices = [sc.Helix(max_offset=32) for _ in range(3)]
    design = sc.Design(helices=helices, grid=sc.square)
    helices[2].roll = 30

    design.draw_strand(0, 0).extension_5p(5, display_length=2.5, display_angle=45) \
        .with_domain_name("ext_5p 1") \
        .move(width).with_domain_name('domain 1') \
        .cross(1).move(-width) \
        .loopout(2, 3).with_domain_name('loopout 1') \
        .move(width) \
        .extension_3p(7).with_domain_name("ext_3p 1") \
        .with_name('strand1') \
        .with_sequence('TATTT' + 'G' * width + 'C' * width + 'T' * 3 + 'G' * width + 'ATAAAAA') \
        .with_modification_5p(sc.cy5_5p) \
        .with_modification_3p(sc.cy3_3p)

    design.draw_strand(0, 24).extension_5p(5, display_length=3.5, display_angle=60) \
        .with_domain_name("ext_5p 2") \
        .move(-width).cross(1).move(width).loopout(2, 3).move(-width) \
        .extension_3p(7).with_domain_name("ext_3p 2") \
        .with_name('strand2') \
        .with_sequence('ATAAA' + 'C' * width + 'G' * width + 'A' * 3 + 'C' * width + 'TATTTTT') \

    return design


if __name__ == '__main__':
    main()
