import scadnano as sc

"""
{
  "color": {"r": 247, "g": 147, "b": 30},
  "dna_sequence": "CCACCAGCAGTACCGAACGAAAGCCCTAATTAACACCGCCTG",
  "substrands": [
    {"helix_idx": 2, "right": false, "start": 289, "end": 299},
    {"helix_idx": 3, "right": true, "start": 289, "end": 310},
    {"helix_idx": 2, "right": false, "start": 299, "end": 310}
  ]
},

{
  "color": {"r": 247, "g": 147, "b": 30},
  "dna_sequence": "???????CCCTAAAACATCGC???????CATTAAAAATACCG",
  "substrands": [
    {"helix_idx": 3, "right": true, "start": 305, "end": 326},
    {"helix_idx": 2, "right": false, "start": 305, "end": 326}
  ]
},
"""

"""
helices = [sc.Helix(0, 330), sc.Helix(1, 330), sc.Helix(2, 330), sc.Helix(3, 330)]

    s1_left_ss0 = sc.Substrand(2, sc.left, 289, 299)
    s1_ss1 = sc.Substrand(3, sc.right, 289, 310)
    s1_right_ss0 = sc.Substrand(2, sc.left, 299, 310)
    s1 = sc.Strand([s1_left_ss0, s1_ss1, s1_right_ss0])

    s2_ss1 = sc.Substrand(3, sc.right, 305, 326)
    s2_ss0 = sc.Substrand(2, sc.left, 305, 326)
    s2 = sc.Strand([s2_ss1, s2_ss0])

    strands = [s1, s2]
    design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)

    return design
"""

def main():
    helices = [sc.Helix(0, 25), sc.Helix(1, 25)]

    s1_left_ss0 = sc.Domain(0, sc.reverse, 0, 5)
    s1_ss1 = sc.Domain(0, sc.forward, 0, 15)
    s1_right_ss0 = sc.Domain(0, sc.reverse, 5, 15)
    s1 = sc.Strand([s1_left_ss0, s1_ss1, s1_right_ss0])

    s2_ss1 = sc.Domain(0, sc.forward, 10, 20)
    s2_ss0 = sc.Domain(0, sc.reverse, 10, 20)
    s2 = sc.Strand([s2_ss1, s2_ss0])

    strands = [s1, s2]
    design = sc.DNADesign(helices=helices, strands=strands, grid=sc.square)

    return design

if not sc.in_browser() and __name__ == '__main__':
    design = main()
    design.write_scadnano_file(directory='output_designs')
