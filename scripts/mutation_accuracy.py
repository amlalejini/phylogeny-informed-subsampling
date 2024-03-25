import csv
import argparse
import os
import statistics


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
    while True:
        evaluated = False
        if len(phylogeny_dict[current_id]['training_cases_evaluated']) > 0:
            evaluated = phylogeny_dict[current_id]['training_cases_evaluated'][training_case_id]
        if evaluated:
            estimate.SetEstimate(
                score = phylogeny_dict[current_id]['phenotype'][training_case_id],
                source = current_id,
                taxon_dist = dist,
                mut_dist = None # TODO - Write function to compute mutational distance
            )
            break
        else:
            current_id = phylogeny_dict[current_id]['ancestor']
            dist += 1
            if current_id in root_ids:
                break

    return estimate

# First comparison
def exhaustive_would_be_estimation(phylogeny_dict, root_ids, extant_ids):
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
            estimates[extant_id].append(estimate_score_ancestor_based(extant_id, i, phylogeny_dict, root_ids))

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
                "mut_dist": None,
                "taxon_dist": distance
            }
            ancestor_accuracies.append(ancestor_accuracy_info)

            if current_id in root_ids:
                break
            distance += 1
            # mut_dist += TODO - come back when we have mutation distance implemented
            current_id = current_taxon_info['ancestor']

        accuracies[extant_id] = ancestor_accuracies

    return accuracies

# Third comparison (WIP)
def accuracy_vs_mut_distance(phylogeny_dict, root_ids, extant_ids):
    """
    For each extant id, collect accuracy and mutation distance over time for all training cases
    {
        "extant_id": [{"accuracy": [...per-training case...], "mut_dist": ...}, ...],
        ...
    }
    """
    accuracy_mut_distance = {}

    for extant_id in extant_ids:
        ancestor_accuracies = [] # List of dictionaries
        # -- Localize extant Info --
        extant_taxon_info = phylogeny_dict[extant_id]
        extant_scores = extant_taxon_info['training_cases_true_scores']
        num_training_cases = len(extant_scores)

        # -- Iteratively walk up ancestors in phylogeny and track time --
        current_id = extant_id
        distance = 0 # Taxon distance
        mut_dist = 0 # Mutation distance
        time = 0 # Time
        while True:
            # Localize current taxon info
            current_taxon_info = phylogeny_dict[current_id]
            current_scores = current_taxon_info["training_cases_true_scores"]
            # Compute accuracy / distance information
            ancestor_accuracy_info = {
                "abs_score_diff": [abs(extant_scores[i] - current_scores[i]) for i in range(num_training_cases)],
                "mut_dist": mut_dist,
            }
            ancestor_accuracies.append(ancestor_accuracy_info)

            if current_id in root_ids:
                break
            distance += 1
            # mut_dist += TODO - come back when we have mutation distance implemented
            current_id = current_taxon_info['ancestor']

        accuracy_mut_distance[extant_id] = ancestor_accuracies

    return accuracy_mut_distance

def parse_list(list_string):
    list_string = list_string.strip("[]").strip()
    if list_string == '':
        return []
    return list_string.split(',')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate mutation accuracy statistics")
    parser.add_argument("--phylo", type=str,  help="Path to phylogeny to analyze")

    args = parser.parse_args()
    phylo_file = args.phylo

    data = read_csv(phylo_file)

    snapshot_update = os.path.basename(phylo_file).split(".")[0]
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
    
    would_be_estimations = exhaustive_would_be_estimation(phylogeny_dict, roots, extant_ids)
    # ---- Output would-be estimations ----
    # (1) Rows = per-extant taxon:
    #     update, summary statistics on accuracy, extant_id
    output_info = []
    # per_taxon_estimates = {}
    for extant_id in would_be_estimations:
        extant_true_scores = phylogeny_dict[extant_id]["training_cases_true_scores"]
        extant_estimations = would_be_estimations[extant_id]
        update = snapshot_update
        num_successes = sum(int(est.estimate_success) for est in extant_estimations)
        # Get distribution of estimation errors
        est_errors = [abs(extant_true_scores[i] - extant_estimations[i].estimated_score) for i in range(len(extant_true_scores)) if extant_estimations[i].estimate_success]
        # Get distribution of estimation taxon distances
        est_dists = [extant_estimations[i].taxon_distance for i in range(len(extant_true_scores))]
        row_info = {
            "update": snapshot_update,
            "num_est_successes": num_successes,
            "extant_id": extant_id,
            "error_mean": statistics.mean(est_errors),
            "error_median": statistics.median(est_errors),
            "error_variance": statistics.variance(est_errors),
            "est_taxon_dist_mean": statistics.mean(est_dists),
            "est_taxon_dist_median": statistics.median(est_dists),
            "est_taxon_dist_variance": statistics.variance(est_dists)
        }
        output_info.append(row_info)
        
    write_csv("taxon_est_info.csv", output_info)
    
    ancestor_accuracy_results = ancestor_vs_extant_scores(phylogeny_dict, roots, extant_ids)
    # ---- Output ancestor vs. extant scores ----
    # (1) Rows = per-extant taxon:
    #     update, summary statistics on accuracy, extant_id
    ancestor_output_info = []
    for extant_id in ancestor_accuracy_results:
        ancestor_accuracies = ancestor_accuracy_results[extant_id]
        extant_true_scores = phylogeny_dict[extant_id]["training_cases_true_scores"]
        update = snapshot_update
        num_cases = len(extant_true_scores)
        for training_case, accuracy_info in enumerate(ancestor_accuracies):
            abs_score_diff = accuracy_info["abs_score_diff"]
            taxon_distance = accuracy_info["taxon_dist"]
            row_info = {
                "update": snapshot_update,
                "extant_id": extant_id,
                "training_case": training_case,
                "accuracy_mean": statistics.mean(abs_score_diff),
                "accuracy_median": statistics.median(abs_score_diff),
                "accuracy_variance": statistics.variance(abs_score_diff),
                "taxon_distance": taxon_distance
            }
            ancestor_output_info.append(row_info)

    write_csv("ancestor_vs_extant_score_info.csv", ancestor_output_info)
