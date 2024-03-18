import csv

class SearchResult:
    def __init__(
        self,
        success = False,
        score = None,
        source = None,
        taxon_dist = None,
        mut_dist = None,
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

def exhaustive_would_be_estimation(phylogeny_dict, root_ids, extant_ids):
    """
    Given phylogeny, roots, and a set of extant ids,
    exhaustively estimate every training case for every extant id.
    """
    # accuracies = []
    estimates = []
    # Loop over extant taxa
    for extant_id in extant_ids:
        taxon_info = phylogeny_dict[extant_id]
        num_training_cases = len(taxon_info['training_cases_true_scores'])
        for i in range(num_training_cases):
            estimates.append(estimate_score_ancestor_based(extant_id, i, phylogeny_dict, root_ids))

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
        ancestor_accuracies = []
        # extant_taxon_info = phylogeny_dict[extant_id]
        # extant_scores = extant_taxon_info['training_cases_true_scores']
        # num_training_cases = len(extant_scores)

        current_id = extant_id
        distance = 0
        num_training_cases = 0  # Define the variable "num_training_cases"
        while True:
            current_taxon_info = phylogeny_dict[current_id]
            # --- Bookmark ---

            ancestor_accuracies.extend(current_taxon_info['phenotype'])
            distance.extend([distance] * num_training_cases)


            if current_taxon_info['ancestor'] is None:
                break
            distance += 1
            current_id = current_taxon_info['ancestor']

    return accuracies


def parse_list(list_string):
    list_string = list_string.strip("[]").strip()
    if list_string == '':
        return []
    return list_string.split(',')

if __name__ == "__main__":
    input_file_path = "/home/sansonm/research_ws/phylogeny-informed-subsampling/output/phylo_1000.csv"

    data = read_csv(input_file_path)

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
    print(f"Would-Be Estimation Values: {would_be_estimations}")
    # proportion_would_be_estimation = sum(would_be_estimations) / len(would_be_estimations)
    # print(f"Would-Be Estimation: {proportion_would_be_estimation:.5f}")
    print()

    ancestor_accuracy = ancestor_vs_extant_scores(phylogeny_dict, roots, extant_ids)
    print(f"Ancestor to Extant Score Accuracy Values: {ancestor_accuracy}")
    # percent_ancestor_accuracy = sum(ancestor_accuracy) / len(ancestor_accuracy)
    # print(f"Ancestor to Extant Score Accuracy: {percent_ancestor_accuracy:.5f}")
