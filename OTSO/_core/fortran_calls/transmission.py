import multiprocessing as mp
import json

from ..utils import mhd_utils
from ..utils import fortran_data_utils
from ..data_classes.transmission_data import TransmissionData
from ..utils import transmission_utils as tu
from ..utils import cpus_utils as cpu_util


def FortranTransmission(Data: list, TransmissionDataInstance: TransmissionData, queue: mp.Queue, cpus, JsonFile, lock) -> None:

    cpu_util.set_process_affinity(cpus)

    if TransmissionDataInstance.model[1] == 99:
        mhd_utils.MHDinitialise(TransmissionDataInstance.MHDfile)

    for x in Data:

        Station = x[0]

        FortranData = fortran_data_utils.prepare_fortran_transmission(
            x, TransmissionDataInstance
        )

        R_high = FortranData.startrigidity
        R_low = FortranData.endrigidity

        Transmissiondf = tu.transmission_full(
                FortranData,
                R_high,
                R_low,
                Station,
                TransmissionDataInstance,
            )

        record = { 
              "station": Station,
              "transmission": Transmissiondf.to_dict(orient="split"),
              }

        with lock:
          with open(JsonFile, "a", encoding="utf-8") as f:
              json.dump(record, f)
              f.write("\n")

        queue.put(1)

    return


