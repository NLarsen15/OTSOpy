from datetime import date

def READMECoordtrans(CoordIN, CoordOUT, Printtime) -> str:

    result = []

    today = date.today()
    result.append(f"\n")
    result.append(f"Date of OTSO computation: {today}\n")
    result.append(f"Total computation time: {Printtime} seconds\n\n")
    result.append(f"Input Coordinate System:\n{CoordIN}\n\n")
    result.append(f"Output Coordinate System:\n{CoordOUT}\n\n")
    

    return "".join(result)