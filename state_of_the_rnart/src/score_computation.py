import os


METRICS = "RMSD,P-VALUE,BARNABA,TM-score,GDT-TS,INF,CAD,lDDT,MCQ"

COMMAND = (
    "docker run -it -v ${PWD}/data/:/app/data "
    "-v ${PWD}/tmp:/tmp sayby77/rnadvisor --pred_path $PRED_PATH "
    "--native_path $NATIVE_PATH --result_path $RESULT_PATH "
    f"--all_scores={METRICS}"
)


class ScoreComputation:
    def __init__(self, pred_dir: str, native_dir: str, out_dir: str):
        """
        Initialise the different parameters to run RNAdvisor
        :param pred_dir: directory where are stored the different predictions
        :param native_dir: directory where are stored the native pdb files
        :param out_dir: where to save the output
        """
        self.pred_dir = pred_dir
        self.native_dir = native_dir
        self.out_dir = out_dir

    def run(self):
        rnas = [name for name in os.listdir(self.native_dir) if name.endswith(".pdb")]
        for rna in rnas:
            in_native = os.path.join(self.native_dir, rna)
            # Prediction path should be a directory and not a file
            in_pred = os.path.join(self.pred_dir, rna.replace(".pdb", ""))
            # Save to a csv format
            out_pred = os.path.join(self.out_dir, rna.replace(".pdb", ".csv"))
            self.run_rnadvisor(in_native, in_pred, out_pred)

    def run_rnadvisor(self, in_native: str, in_pred: str, out_pred: str):
        """
        Run the RNAdvisor tool to get the metrics and evaluate the structures
        :param in_native: path to a pdb file with the native structure
        :param in_pred: path to a directory with the different predictions
        :param out_pred: path to save the computed metrics
        """
        command = (
            COMMAND.replace("$PRED_PATH", in_pred)
            .replace("$NATIVE_PATH", in_native)
            .replace("$RESULT_PATH", out_pred)
        )
        os.system(command)


if __name__ == "__main__":
    params = {
        "pred_dir": os.path.join("data", "input", "PREDS"),
        "native_dir": os.path.join("data", "input", "NATIVE"),
        "out_dir": os.path.join("data", "output"),
    }
    score_computation = ScoreComputation(**params)
    score_computation.run()
