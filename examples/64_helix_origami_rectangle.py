import origami_rectangle as rect
import scadnano as sc


def create_design() -> sc.Design:
    return rect.create(num_helices=64, num_cols=6)


if __name__ == '__main__':
    design = create_design()
    design.write_scadnano_file(directory='output_designs')
