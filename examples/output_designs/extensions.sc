{
  "version": "0.17.5",
  "grid": "square",
  "helices": [
    {"max_offset": 32, "grid_position": [0, 0]},
    {"max_offset": 32, "grid_position": [0, 1]},
    {"max_offset": 32, "grid_position": [0, 2], "roll": 30}
  ],
  "modifications_in_design": {
    "/3Cy3Sp/": {
      "display_text": "Cy3",
      "idt_text": "/3Cy3Sp/",
      "display_connector": false,
      "location": "3'"
    },
    "/5Cy5/": {
      "display_text": "Cy5",
      "idt_text": "/5Cy5/",
      "display_connector": false,
      "location": "5'"
    }
  },
  "strands": [
    {
      "name": "strand1",
      "color": "#f74308",
      "sequence": "TATTTGGGGGGGGCCCCCCCCTTTGGGGGGGGATAAAAA",
      "domains": [
        {"extension_num_bases": 5, "display_length": 2.5, "name": "ext_5p 1"},
        {"name": "domain 1", "helix": 0, "forward": true, "start": 0, "end": 8},
        {"helix": 1, "forward": false, "start": 0, "end": 8},
        {"loopout": 3, "name": "loopout 1"},
        {"helix": 2, "forward": true, "start": 0, "end": 8},
        {"extension_num_bases": 7, "name": "ext_3p 1"}
      ],
      "5prime_modification": "/5Cy5/",
      "3prime_modification": "/3Cy3Sp/"
    },
    {
      "name": "strand2",
      "color": "#57bb00",
      "sequence": "ATAAACCCCCCCCGGGGGGGGAAACCCCCCCCTATTTTT",
      "domains": [
        {"extension_num_bases": 5, "display_length": 3.5, "display_angle": 60, "name": "ext_5p 2"},
        {"helix": 0, "forward": false, "start": 16, "end": 24},
        {"helix": 1, "forward": true, "start": 16, "end": 24},
        {"loopout": 3},
        {"helix": 2, "forward": false, "start": 16, "end": 24},
        {"extension_num_bases": 7, "name": "ext_3p 2"}
      ]
    }
  ]
}