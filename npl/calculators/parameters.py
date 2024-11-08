# This file contains the parameters for the calculators
# The parameters are stored in a dictionary with the stoichiometry as the key
# The value is a dictionary with the following keys:
# - symbols: a list of the symbols in the stoichiometry
# - data: a dictionary with the parameters for each symbol
# - reference: the reference for the parameters

top_parameters = {
    "Pt70Au70": {
        "AuPt": 19,
        "AuAu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Au"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Au": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -619,
                "cn=7": -377,
                "cn=8": -256,
                "cn=9": -256,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pd70Au70": {
        "AuPd": -13,
        "AuAu": 0,
        "PdPd": 0,
        "symbols": ["Pd", "Au"],
        "data": {
            "Pd": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Au": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -404,
                "cn=7": -301,
                "cn=8": -200,
                "cn=9": -200,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt70Ag70": {
        "AgPt": 11,
        "AgAg": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Ag"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Ag": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -625,
                "cn=7": -336,
                "cn=8": -195,
                "cn=9": -195,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pd70Ag70": {
        "AgPd": -1,
        "AgAg": 0,
        "PdPd": 0,
        "symbols": ["Pd", "Ag"],
        "data": {
            "Pd": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Ag": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -361,
                "cn=7": -289,
                "cn=8": -163,
                "cn=9": -163,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt70Cu70": {
        "CuPt": -35,
        "CuCu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Cu"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Cu": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -27,
                "cn=7": 182,
                "cn=8": 344,
                "cn=9": 344,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pd70Cu70": {
        "CuPd": -26,
        "CuCu": 0,
        "PdPd": 0,
        "symbols": ["Pd", "Cu"],
        "data": {
            "Pd": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Cu": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 95,
                "cn=7": 147,
                "cn=8": 183,
                "cn=9": 183,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt151Au50": {
        "AuPt": 21,
        "AuAu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Au"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Au": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -507,
                "cn=7": -543,
                "cn=8": -431,
                "cn=9": -431,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt101Au100": {
        "AuPt": 21,
        "AuAu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Au"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Au": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -530,
                "cn=7": -492,
                "cn=8": -335,
                "cn=9": -335,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt51Au150": {
        "AuPt": 15,
        "AuAu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Au"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Au": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -558,
                "cn=7": -547,
                "cn=8": -259,
                "cn=9": -259,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt151Ag50": {
        "AgPt": 32,
        "AgAg": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Ag"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Ag": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -396,
                "cn=7": -380,
                "cn=8": -237,
                "cn=9": -237,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt101Ag100": {
        "AgPt": 16,
        "AgAg": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Ag"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Ag": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -499,
                "cn=7": -466,
                "cn=8": -308,
                "cn=9": -308,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt51Ag150": {
        "AgPt": 7,
        "AgAg": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Ag"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Ag": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": -408,
                "cn=7": -511,
                "cn=8": -240,
                "cn=9": -240,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt151Cu50": {
        "CuPt": -25,
        "CuCu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Cu"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Cu": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 267,
                "cn=7": 342,
                "cn=8": 372,
                "cn=9": 372,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt101Cu100": {
        "CuPt": -43,
        "CuCu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Cu"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Cu": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 15,
                "cn=7": 208,
                "cn=8": 325,
                "cn=9": 325,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    },
    "Pt51Cu150": {
        "CuPt": -54,
        "CuCu": 0,
        "PtPt": 0,
        "symbols": ["Pt", "Cu"],
        "data": {
            "Pt": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 0,
                "cn=7": 0,
                "cn=8": 0,
                "cn=9": 0,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            },
            "Cu": {
                "cn=0": 0,
                "cn=1": 0,
                "cn=2": 0,
                "cn=3": 0,
                "cn=4": 0,
                "cn=5": 0,
                "cn=6": 132,
                "cn=7": 184,
                "cn=8": 259,
                "cn=9": 259,
                "cn=10": 0,
                "cn=11": 0,
                "cn=12": 0
            }
        },
        "reference": "L. Vega Mater. Adv., 2021, 2, 6589-6602"
    }
}