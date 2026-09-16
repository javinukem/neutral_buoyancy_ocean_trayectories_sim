
def file_name(ellpise:bool, square:bool,intx:int, inty:int, Rx:float, Ry:float, Nabv: int, kAmax, PRES) -> str:
    if ellipse == 1:
        pathname = 'E_IX{}IY{}_RX{}RY{}_NABV{}_kA{}P0{}'.format(intx,inty,Rx,Ry,Nabv,kAmax,PRES[kAmax,0,0])
    if square == 1:
        pathname = 'S_IX{}IY{}_RX{}RY{}_NABV{}_kA{}P0{}'.format(intx,inty,Rx,Ry,Nabv,kAmax,PRES[kAmax,0,0])
    return pathname
