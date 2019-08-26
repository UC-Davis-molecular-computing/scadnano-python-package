import origami_rectangle as rect
import scadnano as sc


def main():
    design = rect.create(num_helices=16, num_cols=28, seam_left_column=11, assign_seq=False,
                         num_flanking_helices=1)
    return design


if not sc.in_browser() and __name__ == '__main__':
    design = main()
    design.write_file(directory='output_designs')
