{
  "version": "0.19.0",
  "grid": "square",
  "helices": [
    {"grid_position": [0, 0]},
    {"grid_position": [0, 1]},
    {"grid_position": [0, 2]},
    {"grid_position": [0, 3]},
    {"grid_position": [0, 4]},
    {"grid_position": [0, 5]},
    {"grid_position": [0, 6]},
    {"grid_position": [0, 7]}
  ],
  "modifications_in_design": {
    "/5Biosg/": {
      "display_text": "B",
      "vendor_code": "/5Biosg/",
      "display_connector": false,
      "location": "5'"
    },
    "/iBiodT/": {
      "display_text": "B",
      "vendor_code": "/iBiodT/",
      "display_connector": false,
      "location": "internal",
      "allowed_bases": ["T"]
    },
    "/3Cy3Sp/": {
      "display_text": "Cy3",
      "vendor_code": "/3Cy3Sp/",
      "display_connector": false,
      "location": "3'"
    },
    "/iCy3/": {
      "display_text": "Cy3",
      "vendor_code": "/iCy3/",
      "display_connector": false,
      "location": "internal"
    }
  },
  "strands": [
    {
      "color": "#f74308",
      "sequence": "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
      "domains": [
        {"helix": 0, "forward": true, "start": 0, "end": 16, "deletions": [11, 12]},
        {"helix": 1, "forward": false, "start": 0, "end": 16, "deletions": [12], "insertions": [[4, 1]]},
        {"helix": 2, "forward": true, "start": 0, "end": 16},
        {"helix": 3, "forward": false, "start": 0, "end": 16},
        {"helix": 4, "forward": true, "start": 0, "end": 16},
        {"helix": 5, "forward": false, "start": 0, "end": 16},
        {"helix": 6, "forward": true, "start": 0, "end": 16},
        {"helix": 7, "forward": false, "start": 0, "end": 16}
      ],
      "5prime_modification": "/5Biosg/",
      "3prime_modification": "/3Cy3Sp/",
      "internal_modifications": {"5": "/iCy3/", "10": "/iBiodT/", "21": "/iCy3/", "26": "/iBiodT/", "37": "/iCy3/", "42": "/iBiodT/", "53": "/iCy3/", "58": "/iBiodT/", "69": "/iCy3/", "74": "/iBiodT/", "85": "/iCy3/", "90": "/iBiodT/", "101": "/iCy3/", "106": "/iBiodT/", "117": "/iCy3/", "122": "/iBiodT/"}
    }
  ]
}