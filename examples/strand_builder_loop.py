import scadnano as sc

def create_design() -> sc.Design:
    num_helices = 32
    helices = [sc.Helix(max_offset=200) for _ in range(num_helices)]
    design = sc.Design(helices=helices, grid=sc.square)
    strand_builder = design.draw_strand(0, 0)
    for helix in range(num_helices):
        # move forward if on an even helix, otherwise move in reverse
        move_distance = 200 if helix % 2 == 0 else -200
        strand_builder.move(move_distance)
        if helix < 31: # crossover to next helix, unless it's the last helix
            strand_builder.cross(helix + 1)
    strand_builder.as_scaffold()
    return design


if __name__ == '__main__':
    design = create_design()
    design.write_scadnano_file(directory='output_designs')