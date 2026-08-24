import json
from pathlib import Path

def get_citations(citation_keys):
    """
    Extract citations from a JSON citation database.

    Parameters:
        json_file (str): Path to JSON citation file
        citation_keys (list): List of citation IDs to retrieve

    Returns:
        list: List of citation dictionaries
    """

    json_file = Path(__file__).parent / "citations.json"

    with open(json_file, "r") as f:
        citation_db = json.load(f)

    citations = []

    for key in citation_keys:
        if key in citation_db["citations"]:
            citations.append(citation_db["citations"][key])
        else:
            print(f"Warning: Citation '{key}' not found in database.")

    return citations


def format_citations(citations):
    """Convert citations into a printable string."""
    
    output = []
    output.append("References Used")
    output.append("================")

    for citation in citations:
        output.append("")
        output.append(citation["citation"])
        output.append(f"DOI: {citation['doi']}")

    return "\n".join(output)


def generate_citation_array(Data):

    Data.citationlist.append("Larsen2023")

    get_internalcitations(Data)
    get_externalcitations(Data)
    get_integrationcitations(Data)
    get_magnetopausecitations(Data)
    get_bobergcitation(Data)

    citations = get_citations(Data.citationlist)

    Data.citationstring = format_citations(citations)

def generate_citation_mag_array(Data):

    Data.citationlist.append("Larsen2023")

    get_internalcitations(Data)
    get_externalcitations(Data)
    get_bobergcitation(Data)

    citations = get_citations(Data.citationlist)

    Data.citationstring = format_citations(citations)


def get_internalcitations(Data):

    if Data.internalmag == "IGRF":
        Data.citationlist.append("IGRF14")
    if Data.internalmag == "CHAOS":
        Data.citationlist.append("CHAOS")


def get_externalcitations(Data):

    if Data.externalmag == "TSY87short" or Data.externalmag == "TSY87long":
        Data.citationlist.append("T87")

    elif Data.externalmag == "TSY89a":
        Data.citationlist.append("T89")

    elif Data.externalmag == "TSY89c":
        Data.citationlist.append("T89")
        Data.citationlist.append("T89c")

    elif Data.externalmag == "TSY89_refit":
        Data.citationlist.append("T89")
        Data.citationlist.append("T89_refit")

    elif Data.externalmag == "T96":
        Data.citationlist.append("T96")

    elif Data.externalmag == "T01":
        Data.citationlist.append("T01a")
        Data.citationlist.append("T01b")

    elif Data.externalmag == "T01S":
        Data.citationlist.append("T01S")

    elif Data.externalmag == "T04":
        Data.citationlist.append("T04")

    elif Data.externalmag == "T15":
        Data.citationlist.append("T15")

    elif Data.externalmag == "T16":
        Data.citationlist.append("T16")

def get_integrationcitations(Data):

    if Data.intmodel == "4RK":
        Data.citationlist.append("RKa")
        Data.citationlist.append("RKb")

    elif Data.intmodel == "Boris":
        Data.citationlist.append("Boris1970")

    elif Data.intmodel == "Boris-Buneman":
        Data.citationlist.append("Boris1970")

    elif Data.intmodel == "HC":
        Data.citationlist.append("HigueraCary2017")

    elif Data.intmodel == "Vay":
        Data.citationlist.append("Vay2008")

def get_magnetopausecitations(Data):

    if Data.magnetopause == "Kobel":
        Data.citationlist.append("Kobel1994")

    elif Data.magnetopause == "aFormisano":
        Data.citationlist.append("Formisano79")

    elif Data.magnetopause == "Sibeck":
        Data.citationlist.append("Sibeck1991")

    elif Data.magnetopause == "Lin":
        Data.citationlist.append("Lin2010")

def get_bobergcitation(Data):

    if Data.boberg:
        Data.citationlist.append("Boberg1995")