import csv
import argparse
import os
import statistics
from collections import deque
import utilities as utils

run_cfg_fields = [
    "SEED",
    "POP_SIZE",
    "STOP_MODE",
    "MAX_GENS",
    "MAX_EVALS",
    "POP_INIT_MODE",
    "SELECTION",
    "TOURNAMENT_SIZE",
    "PROBLEM",
    "TESTING_SET_PATH",
    "TRAINING_SET_PATH",
    "ANCESTOR_FILE_PATH",
    "EVAL_MODE",
    "EVAL_FIT_EST_MODE",
    "EVAL_MAX_PHYLO_SEARCH_DEPTH",
    "EVAL_CPU_CYCLES_PER_TEST",
    "NUM_COHORTS",
    "TEST_DOWNSAMPLE_RATE",
    "MAX_ACTIVE_THREAD_CNT",
    "MAX_THREAD_CAPACITY",
    "PRG_MIN_FUNC_CNT",
    "PRG_MAX_FUNC_CNT",
    "PRG_MIN_FUNC_INST_CNT",
    "PRG_MAX_FUNC_INST_CNT",
    "PRG_INST_MIN_ARG_VAL",
    "PRG_INST_MAX_ARG_VAL",
    "MUT_RATE_INST_ARG_SUB",
    "MUT_RATE_INST_SUB",
    "MUT_RATE_INST_INS",
    "MUT_RATE_INST_DEL",
    "MUT_RATE_SEQ_SLIP",
    "MUT_RATE_FUNC_DUP",
    "MUT_RATE_FUNC_DEL",
    "MUT_RATE_INST_TAG_BF",
    "MUT_RATE_FUNC_TAG_BF",
    "MUT_RATE_INST_TAG_SINGLE_BF",
    "MUT_RATE_FUNC_TAG_SINGLE_BF",
    "MUT_RATE_INST_TAG_SEQ_RAND",
    "MUT_RATE_FUNC_TAG_SEQ_RAND"
]

per_training_case_cfg = {
    "SEED"
}

def write_csv(output_path, summary_dict):
    # (1) What's the header?
    header = list(summary_dict[0].keys())
    header.sort()
    # (2) Collect all lines as strings
    lines = []
    for info in summary_dict:
        line_header_info = sorted(list(info.keys()))
        if line_header_info != header:
            print("Header mismatch!")
            exit(-1)
        line = ",".join([str(info[field]) for field in header])
        lines.append(line)
    out_content = ",".join(header) + "\n"
    out_content += "\n".join(lines)

    with open(output_path, "w") as fp:
        fp.write(out_content)

    return header

def append_csv(output_path, out_lines, field_order):
    lines = []
    for info in out_lines:
        line = ",".join([str(info[field]) for field in field_order])
        lines.append(line)
    out_content = "\n" + "\n".join(lines)
    with open(output_path, "a") as fp:
        fp.write(out_content)

class SearchResult:
    def __init__(
        self,
        success = False,
        score = None,
        source = None,
        taxon_dist = None,
        mut_dist = None
    ):
        self.estimate_success = success # Estimation successful?
        self.estimated_score = score   # Score from source of estimate
        self.estimate_source = source   # Source id for the estimate
        self.taxon_distance = taxon_dist    # How many steps in the phylogeny away is the source of the estimate?
        self.mutation_distance = mut_dist # How mutationally distant is the source of the estimate?

    def SetEstimate(self, score, source, taxon_dist, mut_dist):
        self.estimate_success = True
        self.estimated_score = score
        self.estimate_source = source
        self.taxon_distance = taxon_dist
        self.mutation_distance = mut_dist

    def SetFail(self):
        self.estimate_success = False
        self.estimated_score = None
        self.estimate_source = None
        self.taxon_distance = None
        self.mutation_distance = None

    def __str__(self) -> str:
        return f"Success:{self.estimate_success}, Score:{self.estimated_score}, Source:{self.estimate_source}, Taxon dist:{self.taxon_distance}, Mut dist:{self.mutation_distance}"

def read_csv(file_path):
    data = []
    with open(file_path, "r", newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            data.append(row)
    return data

def estimate_score_ancestor_based(taxon_id, training_case_id, phylogeny_dict, root_ids):
    """
    Given a taxon_id, training_case_id, the phylogeny, and roots, return SearchResult
    using ancestor-based estimation.
    """
    # taxon_info = phylogeny_dict[taxon_id]
    estimate = SearchResult()

    # Ancestor-based estimation
    current_id = taxon_id
    dist = 0
    mut_dist = 0
    while True:
        evaluated = False
        if len(phylogeny_dict[current_id]['training_cases_evaluated']) > 0:
            evaluated = phylogeny_dict[current_id]['training_cases_evaluated'][training_case_id]
        if evaluated:
            estimate.SetEstimate(
                score = phylogeny_dict[current_id]['phenotype'][training_case_id],
                source = current_id,
                taxon_dist = dist,
                mut_dist = mut_dist
            )
            break
        else:
            current_id = phylogeny_dict[current_id]['ancestor']
            dist += 1
            mut_dist += phylogeny_dict[current_id]["mutational_dist"]
            if current_id in root_ids:
                break

    return estimate

def estimate_score_relative_based(
    taxon_id,
    training_case_id,
    phylogeny_dict,
    root_ids
):
    # Track discovered taxa ids (starting with current taxon id)
    discovered_taxa = {taxon_id}
    # Use queue to manage breadth first search
    search_queue = deque()
    search_queue.append(
        {"id": taxon_id, "dist": 0, "mut_dist": 0}
    )
    estimate = SearchResult()
    while (len(search_queue) > 0):
        # Grab the current taxon and its distance from search source
        cur_taxon_id = search_queue[0]["id"]
        cur_dist = search_queue[0]["dist"]
        cur_mut_dist = search_queue[0]["mut_dist"]
        # Pop current taxon from the front of the search queue
        search_queue.popleft()
        # Localize relevant taxon info
        taxon_info = phylogeny_dict[cur_taxon_id]
        traits_evaluated = taxon_info["training_cases_evaluated"]
        # Was focal training case evalutaed on this taxon?
        evaluated = len(traits_evaluated) > 0 and traits_evaluated[training_case_id]
        if evaluated:
            estimate.SetEstimate(
                score = taxon_info['phenotype'][training_case_id],
                source = cur_taxon_id,
                taxon_dist = cur_dist,
                mut_dist = cur_mut_dist
            )
            return estimate
        # Training case was not evaluated
        # - Add any ancestors not already searched
        ancestor_taxon_id = taxon_info["ancestor"]
        if ((not cur_taxon_id in root_ids) and (not ancestor_taxon_id in discovered_taxa)):
            discovered_taxa.add(ancestor_taxon_id)
            # If moving _to_ an ancestor, mut distance equals mutation distance stored
            # on current taxon
            search_queue.append(
                {
                    "id": ancestor_taxon_id,
                    "dist": cur_dist + 1,
                    "mut_dist": cur_mut_dist + taxon_info["mutational_dist"]
                }
            )
        # - Add any descendents not already searched
        descendant_taxa_ids = sorted(list(taxon_info["descendents"]))
        for descendant_id in descendant_taxa_ids:
            if not (descendant_id in discovered_taxa):
                # If we're moving to a descendent, mutation distances increases
                # by amount on that descendent
                # Note that we can never go back up from a descendent with an asexual population
                # (so no risk of double counting)
                discovered_taxa.add(descendant_id)
                search_queue.append(
                    {
                        "id": descendant_id,
                        "dist": cur_dist + 1,
                        "mut_dist": cur_mut_dist + phylogeny_dict[descendant_id]["mutational_dist"]
                    }
                )

    return estimate

def exhaustive_would_be_estimation(
    phylogeny_dict,
    root_ids,
    extant_ids,
    estimation_mode = "ancestor"
):
    """
    Given phylogeny, roots, and a set of extant ids,
    exhaustively estimate every training case for every extant id.
    Returns: {
        "extant_id": [list of SearchResults (per training case)]
    }
    """
    # accuracies = []
    estimates = {}
    # Loop over extant taxa
    for extant_id in extant_ids:
        estimates[extant_id] = []
        taxon_info = phylogeny_dict[extant_id]
        num_training_cases = len(taxon_info['training_cases_true_scores'])
        for i in range(num_training_cases):
            if estimation_mode == "ancestor":
                estimates[extant_id].append(estimate_score_ancestor_based(extant_id, i, phylogeny_dict, root_ids))
            elif estimation_mode == "relative":
                estimates[extant_id].append(estimate_score_relative_based(extant_id, i, phylogeny_dict, root_ids))
            else:
                print("Unknown estimation mode:", estimation_mode)
                exit(-1)
    return estimates

# Second comparison
def ancestor_vs_extant_scores(phylogeny_dict, root_ids, extant_ids):
    """
    For each extant id, collect accuracy of each of its ancestors for all training cases
    {
        "extant_id": [{"accuracy": [...per-training case...], "mut_dist":, "taxon_dist":...}, ...],
        ...
    }
    """
    accuracies = {}

    for extant_id in extant_ids:
        ancestor_accuracies = [] # List of dictionaries
        # -- Localize extant Info --
        extant_taxon_info = phylogeny_dict[extant_id]
        extant_scores = extant_taxon_info['training_cases_true_scores']
        num_training_cases = len(extant_scores)

        # -- Iteratively walk up ancestors in phylogeny --
        current_id = extant_id
        distance = 0 # Taxon distance
        mut_dist = 0 # Mutation distance
        while True:
            # Localize current taxon info
            current_taxon_info = phylogeny_dict[current_id]
            current_scores = current_taxon_info["training_cases_true_scores"]
            # Compute accuracy / distance information
            ancestor_accuracy_info = {
                "abs_score_diff": [abs(extant_scores[i] - current_scores[i]) for i in range(num_training_cases)],
                "mut_dist": mut_dist,
                "taxon_dist": distance
            }
            ancestor_accuracies.append(ancestor_accuracy_info)

            if current_id in root_ids:
                break
            distance += 1
            mut_dist += current_taxon_info["mutational_dist"]
            current_id = current_taxon_info['ancestor']

        accuracies[extant_id] = ancestor_accuracies

    return accuracies


def parse_list(list_string):
    list_string = list_string.strip("[]").strip()
    if list_string == '':
        return []
    return list_string.split(',')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mutation accuracy statistics")
    # parser.add_argument("--phylo", type=str,  help="Path to phylogeny to analyze")
    parser.add_argument("--data_dir", type=str, nargs="+", help="Where is the base output directory for each run?")
    parser.add_argument("--dump", type=str, help="Where to dump this?", default=".")
    parser.add_argument("--run_dir_id", type=str, help="String contained in all run directories", default="RUN_")

    # Parse command line arguments
    args = parser.parse_args()
    data_dirs = args.data_dir
    dump_dir = args.dump
    run_identifier = args.run_dir_id

    # Verify that the given data directory exits
    if not all([os.path.exists(data_dir) for data_dir in data_dirs]):
        print("Unable to find data directory.")
        exit(-1)

    utils.mkdir_p(dump_dir)

    # Aggregate run directories
    run_dirs = []
    for data_dir in data_dirs:
        run_dirs += [os.path.join(data_dir, run_dir) for run_dir in os.listdir(data_dir) if run_identifier in run_dir]
    print(f"Found {len(run_dirs)} run directories.")

    total_runs = len(run_dirs)
    cur_run_i = 0

    header_extant_taxon_est_info = None
    header_extant_training_case_est_info = None
    header_ancestor_error = None
    header_ancestor_error_per_training_case = None
    header_cfg_file = None
    for run_path in run_dirs:
        cur_run_i += 1
        print(f"Processing ({cur_run_i}/{total_runs}): {run_path}")
        output_path = os.path.join(run_path, "output")

        cfg_info = {}
        ############################################################
        # Extract configs from run_config.csv file
        cfg_path = os.path.join(output_path, "run_config.csv")
        cfg_data = utils.read_csv(cfg_path)
        cmd_params = {}
        for line in cfg_data:
            param = line["parameter"]
            value = line["value"]
            cmd_params[param] = value
            if param in run_cfg_fields:
                cfg_info[param] = value
        print("Run configuration:", cfg_info)
        cmd_params["instruction_set"] = f'"{cmd_params["instruction_set"]}"'
        if header_cfg_file is None:
            header_cfg_file = write_csv(
                os.path.join(dump_dir, "run_cfgs.csv"),
                [cmd_params]
            )
        else:
            append_csv(
                os.path.join(dump_dir, "run_cfgs.csv"),
                [cmd_params],
                header_cfg_file
            )
        ############################################################

        # TODO - Identify phylogeny file(s) to analyze
        phylo_files = [fname for fname in os.listdir(output_path) if "phylo_" in fname and ".csv" in fname]
        final_phylo_file = max(phylo_files, key = lambda x: int(x.split(".")[0].split("_")[-1]))
        print(final_phylo_file)

        # exit(-1)
        # --- make play nice with looping over directories ---
        data = read_csv(os.path.join(output_path, final_phylo_file))

        snapshot_update = os.path.basename(final_phylo_file).split(".")[0].split("_")[-1]
        phylogeny_dict = {}
        roots = set()
        extant_ids = set()

        for row_i in range(len(data)):
            taxon_id = data[row_i]['id']
            phylogeny_dict[taxon_id] = {key:data[row_i][key] for key in data[row_i]}

            evaluated_cases = parse_list(data[row_i]['training_cases_evaluated'])
            phenotype_scores = parse_list(data[row_i]['phenotype'])
            true_scores = parse_list(data[row_i]['training_cases_true_scores'])
            phylogeny_dict[taxon_id]['phenotype'] = list(map(float, phenotype_scores))
            phylogeny_dict[taxon_id]['training_cases_true_scores'] = list(map(float, true_scores))
            phylogeny_dict[taxon_id]['training_cases_evaluated'] = [bool(int(val)) for val in evaluated_cases]
            phylogeny_dict[taxon_id]['mutational_dist'] = int(phylogeny_dict[taxon_id]['mutational_dist'])

            ancestor = phylogeny_dict[taxon_id]['ancestor_list'].strip('[]')
            phylogeny_dict[taxon_id]['ancestor'] = ancestor
            if ancestor == 'NONE':
                roots.add(taxon_id)

            if phylogeny_dict[taxon_id]['destruction_time'] == "inf":
                extant_ids.add(taxon_id)

            phylogeny_dict[taxon_id]['descendents'] = set()

        for taxon_id in phylogeny_dict:
            if taxon_id in roots:
                continue
            ancestor_id = phylogeny_dict[taxon_id]['ancestor']
            phylogeny_dict[ancestor_id]['descendents'].add(taxon_id)

        # Analyze would-be estimations.
        # I.e., exhaustively estimate every training case for every extant taxon.
        would_be_estimations = {}
        would_be_estimations["ancestor"] = exhaustive_would_be_estimation(
            phylogeny_dict,
            roots,
            extant_ids,
            estimation_mode="ancestor"
        )
        would_be_estimations["relative"] = exhaustive_would_be_estimation(
            phylogeny_dict,
            roots,
            extant_ids,
            estimation_mode="relative"
        )
        assert set(would_be_estimations["relative"].keys()) == set(would_be_estimations["ancestor"].keys())
        assert set(would_be_estimations["relative"].keys()) == extant_ids

        # ---- Output would-be estimations ----
        # (1) Rows = per-extant taxon:
        #     update, summary statistics on accuracy, extant_id
        output_info = []
        for est_mode in would_be_estimations:
            for extant_id in extant_ids:
                # -- Append ancestor-based estimation info --
                extant_true_scores = phylogeny_dict[extant_id]["training_cases_true_scores"]
                extant_estimations = would_be_estimations[est_mode][extant_id]
                update = snapshot_update
                num_successes = sum(int(est.estimate_success) for est in extant_estimations)
                # Get distribution of estimation errors
                est_errors = [abs(extant_true_scores[i] - extant_estimations[i].estimated_score) for i in range(len(extant_true_scores)) if extant_estimations[i].estimate_success]
                # Get distribution of estimation taxon distances
                est_dists = [extant_estimations[i].taxon_distance for i in range(len(extant_true_scores))]
                # Get distribution of estimation mutation distances
                est_mut_dists = [extant_estimations[i].mutation_distance for i in range(len(extant_true_scores))]
                row_info = {
                    "update": snapshot_update,
                    "num_est_successes": num_successes,
                    "estimation_mode": est_mode,
                    "extant_id": extant_id,
                    "error_mean": statistics.mean(est_errors),
                    "error_median": statistics.median(est_errors),
                    "error_variance": statistics.variance(est_errors),
                    "est_taxon_dist_mean": statistics.mean(est_dists),
                    "est_taxon_dist_median": statistics.median(est_dists),
                    "est_taxon_dist_variance": statistics.variance(est_dists),
                    "est_mut_dist_mean": statistics.mean(est_mut_dists),
                    "est_mut_dist_median": statistics.median(est_mut_dists),
                    "est_mut_dist_variance": statistics.variance(est_mut_dists)
                }
                for cfg_field in cfg_info:
                    row_info[cfg_field] = cfg_info[cfg_field]

                output_info.append(row_info)
        if header_extant_taxon_est_info is None:
            header_extant_taxon_est_info = write_csv(
                os.path.join(dump_dir, "extant_taxon_est_info.csv"),
                output_info
            )
        else:
            append_csv(
                os.path.join(dump_dir, "extant_taxon_est_info.csv"),
                output_info,
                header_extant_taxon_est_info
            )

        # (2) Rows = per-training case
        output_info = []
        for est_mode in would_be_estimations:
            for extant_id in extant_ids:
                extant_true_scores = phylogeny_dict[extant_id]["training_cases_true_scores"]
                extant_estimations = would_be_estimations[est_mode][extant_id]
                update = snapshot_update
                for training_case_id in range(len(extant_estimations)):
                    error = abs(
                        extant_true_scores[training_case_id] - extant_estimations[training_case_id].estimated_score
                    )
                    row_info = {
                        "update": update,
                        "training_case_id": training_case_id,
                        "estimation_mode": est_mode,
                        "extant_id": extant_id,
                        "est_success": extant_estimations[training_case_id].estimate_success,
                        "error": error,
                        "est_taxon_dist": extant_estimations[training_case_id].taxon_distance,
                        "est_mut_dist": extant_estimations[training_case_id].mutation_distance
                    }
                    for cfg_field in cfg_info:
                        if not cfg_field in per_training_case_cfg: continue
                        row_info[cfg_field] = cfg_info[cfg_field]
                    output_info.append(row_info)

        # header_ancestor_error
        # header_ancestor_error_per_training_case
        if header_extant_training_case_est_info is None:
            header_extant_training_case_est_info = write_csv(
                os.path.join(dump_dir, "extant_training_case_est_info.csv"),
                output_info
            )
        else:
            append_csv(
                os.path.join(dump_dir, "extant_training_case_est_info.csv"),
                output_info,
                header_extant_training_case_est_info
            )

        # ---- Relationship between ancestor distance and extant scores ----
        # "extant_id": [{"accuracy": [...per-training case...], "mut_dist":, "taxon_dist":...}, ...],
        ancestor_accuracy = ancestor_vs_extant_scores(phylogeny_dict, roots, extant_ids)
        assert extant_ids == set(ancestor_accuracy.keys())
        # (1) 1 row = taxon, ancestor dist pair (training cases summarized)
        output_info = []
        for extant_id in extant_ids:
            accuracy_info = ancestor_accuracy[extant_id]
            update = snapshot_update
            for ancestor_info in accuracy_info:
                ancestor_dist = ancestor_info["taxon_dist"]
                mut_dist = ancestor_info["mut_dist"]
                ancestor_errors = ancestor_info["abs_score_diff"]
                row_info = {
                    "update": snapshot_update,
                    "extant_id": extant_id,
                    "error_mean": statistics.mean(ancestor_errors),
                    "error_median": statistics.median(ancestor_errors),
                    "error_variance": statistics.variance(ancestor_errors),
                    "mutation_dist": mut_dist,
                    "taxon_dist": ancestor_dist
                }
                for cfg_field in cfg_info:
                    row_info[cfg_field] = cfg_info[cfg_field]
                output_info.append(row_info)

                # print(ancestor_info)
                # exit(-1)
        # header_ancestor_error_per_training_case
        if header_ancestor_error is None:
            header_ancestor_error = write_csv(
                os.path.join(dump_dir, "ancestor_error.csv"),
                output_info
            )
        else:
            append_csv(
                os.path.join(dump_dir, "ancestor_error.csv"),
                output_info,
                header_ancestor_error
            )

        # (2) Per-training case per-taxon per-ancestor distance
        output_info = []
        for extant_id in extant_ids:
            accuracy_info = ancestor_accuracy[extant_id]
            update = snapshot_update
            for ancestor_info in accuracy_info:
                ancestor_dist = ancestor_info["taxon_dist"]
                mut_dist = ancestor_info["mut_dist"]
                ancestor_errors = ancestor_info["abs_score_diff"]
                for training_case_id in range(len(ancestor_errors)):
                    row_info = {
                        "update": snapshot_update,
                        "extant_id": extant_id,
                        "training_case_id": training_case_id,
                        "error": ancestor_errors[training_case_id],
                        "mutation_dist": mut_dist,
                        "taxon_dist": ancestor_dist
                    }
                    for cfg_field in cfg_info:
                        if not cfg_field in per_training_case_cfg: continue
                        row_info[cfg_field] = cfg_info[cfg_field]
                    output_info.append(row_info)

        if header_ancestor_error_per_training_case is None:
            header_ancestor_error_per_training_case = write_csv(
                os.path.join(dump_dir, "ancestor_error_per_training_case.csv"),
                output_info
            )
        else:
            append_csv(
                os.path.join(dump_dir, "ancestor_error_per_training_case.csv"),
                output_info,
                header_ancestor_error_per_training_case
            )

        # TODO - Relative distance vs extant scores
        # TODO (?) - Is success on some training cases predictive of success on others?
        # TODO - setup to loop over experiment directory structure, operate over multiple phylogeny snapshots
        #   --> This should probably be a separate analysis (it's a different question + this file is already complicated)


