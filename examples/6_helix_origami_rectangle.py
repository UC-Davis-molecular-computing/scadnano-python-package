import origami_rectangle as rect
import scadnano as sc


def create_design() -> sc.Design:
    design = rect.create(num_helices=6, num_cols=10, nick_pattern=rect.staggered, twist_correction_deletion_spacing=3)
    return design


if __name__ == '__main__':
    d = create_design()
    d.write_scadnano_file(directory='output_designs')
