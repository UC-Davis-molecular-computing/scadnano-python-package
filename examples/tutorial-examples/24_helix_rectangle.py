import scadnano as sc


def main() -> None:
    design = create_design()
    design.write_idt_plate_excel_file(use_default_plates=True)
    design.write_scadnano_file()


def create_design() -> sc.Design:
    design = helices_only()

    add_scaffold_precursors(design)
    add_scaffold_crossovers(design)

    add_staple_precursors(design)
    add_staple_crossovers(design)
    add_staple_nicks(design)

    add_twist_correction_deletions(design)
    design.assign_m13_to_scaffold()

    return design


def helices_only() -> sc.Design:
    helices = [sc.Helix(max_offset=288) for _ in range(24)]
    return sc.Design(helices=helices, grid=sc.square)


def add_scaffold_precursors(design: sc.Design) -> None:
    for helix in range(0, 23, 2):  # scaffold goes forward on even helices
        design.strand(helix, 0).move(288).as_scaffold()
    for helix in range(1, 23, 2):  # scaffold goes reverse on odd helices
        design.strand(helix, 288).move(-288).as_scaffold()
    design.strand(23, 288).move(-144).as_scaffold()
    design.strand(23, 144).move(-144).as_scaffold()


def add_scaffold_crossovers(design: sc.Design) -> None:
    for helix in range(1, 23, 2):  # scaffold interior crossovers
        design.add_full_crossover(helix=helix, helix2=helix + 1, offset=144, forward=False)

    for helix in range(0, 23, 2):  # scaffold edges crossovers
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=0, forward=True)
        design.add_half_crossover(helix=helix, helix2=helix + 1, offset=287, forward=True)


def add_staple_precursors(design: sc.Design) -> None:
    staples = [sc.Strand([sc.Domain(helix=helix, forward=helix % 2 == 1, start=0, end=288)])  # noqa
               for helix in range(24)]
    for staple in staples:
        design.add_strand(staple)


def add_staple_crossovers(design: sc.Design) -> None:
    for helix in range(23):
        start_offset = 16 if helix % 2 == 0 else 32
        for offset in range(start_offset, 288, 32):
            if offset != 144:  # skip crossover near seam
                design.add_full_crossover(helix=helix, helix2=helix + 1, offset=offset,
                                          forward=helix % 2 == 1)


def add_staple_nicks(design: sc.Design) -> None:
    for helix in range(24):
        start_offset = 24 if helix % 2 == 0 else 40
        for offset in range(start_offset, 272, 32):
            design.add_nick(helix, offset, forward=helix % 2 == 1)


def add_twist_correction_deletions(design: sc.Design) -> None:
    for helix in range(24):
        for offset in range(19, 286, 48):
            design.add_deletion(helix, offset)




if __name__ == '__main__':
    main()
