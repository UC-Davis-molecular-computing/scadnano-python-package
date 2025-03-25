{
  "version": "0.19.5",
  "groups": {
    "group 0": {
      "position": {"x": 0, "y": 0, "z": 0},
      "grid": "square"
    },
    "group 1": {
      "position": {"x": 0, "y": 3, "z": 0},
      "grid": "square",
      "geometry": {
        "bases_per_turn": 18
      }
    }
  },
  "helices": [
    {"group": "group 0", "grid_position": [0, 0]},
    {"group": "group 1", "grid_position": [0, 0]}
  ],
  "strands": [
    {
      "color": "#f74308",
      "domains": [
        {"helix": 0, "forward": true, "start": 0, "end": 40}
      ]
    },
    {
      "color": "#57bb00",
      "domains": [
        {"helix": 0, "forward": false, "start": 0, "end": 40}
      ]
    },
    {
      "color": "#888888",
      "domains": [
        {"helix": 1, "forward": true, "start": 0, "end": 40}
      ]
    },
    {
      "color": "#32b86c",
      "domains": [
        {"helix": 1, "forward": false, "start": 0, "end": 40}
      ]
    }
  ]
}