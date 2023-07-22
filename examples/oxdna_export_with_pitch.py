import scadnano as sc


def create_design(pitch: float) -> sc.Design:
    helices = [sc.Helix(max_offset=100) for _ in range(2)]
    design = sc.Design(helices=helices, grid=sc.square)

    design.groups[sc.default_group_name].pitch = pitch

    design.draw_strand(0, 0).to(21)
    design.draw_strand(1, 0).to(21)

    return design


if __name__ == '__main__':
    for pitch in [0, 45]:
        d = create_design(pitch)
        d.write_scadnano_file(filename=f'oxdna_export_with_pitch_{pitch}.sc', directory='output_designs')
        d.write_oxdna_files(filename_no_extension=f'oxdna_export_with_pitch_{pitch}', directory='oxdna')
