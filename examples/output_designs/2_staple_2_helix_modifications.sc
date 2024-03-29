{
  "version": "0.19.0",
  "grid": "square",
  "helices": [
    {"grid_position": [0, 0]},
    {"grid_position": [0, 1]}
  ],
  "modifications_in_design": {
    "/5Biosg/": {
      "display_text": "B",
      "vendor_code": "/5Biosg/",
      "display_connector": false,
      "location": "5'"
    },
    "/3Bio/": {
      "display_text": "B",
      "vendor_code": "/3Bio/",
      "display_connector": false,
      "location": "3'"
    },
    "/iBiodT/": {
      "display_text": "B",
      "vendor_code": "/iBiodT/",
      "display_connector": false,
      "location": "internal",
      "allowed_bases": ["T"]
    },
    "/iCy3/": {
      "display_text": "Cy3",
      "vendor_code": "/iCy3/",
      "display_connector": false,
      "location": "internal"
    },
    "/iCy5/": {
      "display_text": "Cy5",
      "vendor_code": "/iCy5/",
      "display_connector": false,
      "location": "internal"
    },
    "/3Cy3Sp/": {
      "display_text": "Cy3",
      "vendor_code": "/3Cy3Sp/",
      "display_connector": false,
      "location": "3'"
    },
    "/5Cy5/": {
      "display_text": "Cy5",
      "vendor_code": "/5Cy5/",
      "display_connector": false,
      "location": "5'"
    }
  },
  "strands": [
    {
      "color": "#0066cc",
      "sequence": "AACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTAACTA",
      "domains": [
        {"helix": 1, "forward": false, "start": 0, "end": 16, "deletions": [12], "insertions": [[2, 1]]},
        {"helix": 0, "forward": true, "start": 0, "end": 32, "deletions": [11, 12, 24], "insertions": [[29, 1]]},
        {"helix": 1, "forward": false, "start": 16, "end": 32, "deletions": [24]}
      ],
      "5prime_modification": "/5Biosg/",
      "3prime_modification": "/3Cy3Sp/",
      "internal_modifications": {"5": "/iCy5/", "32": "/iCy3/"}
    },
    {
      "color": "#f74308",
      "sequence": "AGTTAGTTAGTTAGTTTTAGTTAGTTAGTT",
      "domains": [
        {"helix": 1, "forward": true, "start": 0, "end": 16, "deletions": [12], "insertions": [[2, 1]]},
        {"helix": 0, "forward": false, "start": 0, "end": 16, "deletions": [11, 12]}
      ],
      "5prime_modification": "/5Biosg/",
      "3prime_modification": "/3Cy3Sp/",
      "internal_modifications": {"9": "/iCy3/", "10": "/iBiodT/", "11": "/iCy3/", "12": "/iCy5/", "4": "/iCy3/", "26": "/iCy5/"}
    },
    {
      "color": "#57bb00",
      "sequence": "TTAGTTAGTTAGTTAGTAGTTAGTTAGTTAG",
      "domains": [
        {"helix": 0, "forward": false, "start": 16, "end": 32, "deletions": [24], "insertions": [[29, 1]]},
        {"helix": 1, "forward": true, "start": 16, "end": 32, "deletions": [24]}
      ],
      "5prime_modification": "/5Cy5/",
      "3prime_modification": "/3Bio/",
      "internal_modifications": {"5": "/iCy3/"}
    }
  ]
}