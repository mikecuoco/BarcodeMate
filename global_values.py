BARCODE="BX"
def change_to_molecule_barcode():
    global BARCODE
    BARCODE="MI"
def is_molecule_barcode():
    if BARCODE=="MI":
        return True
    else:
        return False