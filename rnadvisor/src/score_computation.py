import os


METRICS = "RMSD,P-VALUE,INF,DI,MCQ,TM-SCORE,lDDT,rsRNASP"

COMMAND = (
    "docker run -it -v ${PWD}/docker_data/:/app/docker_data "
    "-v ${PWD}/tmp:/tmp sayby77/rnadvisor --pred_path $PRED_PATH "
    "--native_path $NATIVE_PATH --result_path $RESULT_PATH "
    f"--all_scores={METRICS}"
)
import os


METRICS = "RMSD,P-VALUE,INF,DI,MCQ,TM-SCORE,lDDT,rsRNASP"

COMMAND = (
    "docker run -it -v ${PWD}/docker_data/:/app/docker_data "
    "-v ${PWD}/tmp:/tmp sayby77/rnadvisor --pred_path $PRED_PATH "
    "--native_path $NATIVE_PATH --result_path $RESULT_PATH "
    f"--all_scores={METRICS}"
)


if __name__ == "__main__":
    native_path = os.path.join("../docker_data", "input", "NATIVE", "3_solution_0_rpr.pdb")
    pred_path = os.path.join("../docker_data", "input", "PREDS")
    result_path = os.path.join("../docker_data", "output", "rp03_results_bis.csv")
    command = COMMAND.replace("$NATIVE_PATH", native_path).replace("$PRED_PATH", pred_path).\
        replace("$RESULT_PATH", result_path)
    os.system(command)


if __name__ == "__main__":
    native_path = os.path.join("../docker_data", "input", "NATIVE", "3_solution_0_rpr.pdb")
    pred_path = os.path.join("../docker_data", "input", "PREDS")
    result_path = os.path.join("../docker_data", "output", "rp03_results_bis.csv")
    command = COMMAND.replace("$NATIVE_PATH", native_path).replace("$PRED_PATH", pred_path).\
        replace("$RESULT_PATH", result_path)
    os.system(command)
