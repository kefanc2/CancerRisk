from analyzer import Analyzer

def get_max_by_time(tree, num_timepoints):
    # tree: a list of paths, each path is a list of nodes [num_mutations, score, clone id]
    max_t = max([tree[i][-1][0] for i in range(len(tree))])
    out = [tree[0][0][1]] * num_timepoints
    max_clones = [set()] * num_timepoints
    for path in tree:
        for node in path[::-1]:
            t = min(node[0] * num_timepoints // max_t, num_timepoints-1)
            score = node[1]
            for i in range(int(t), num_timepoints):
                if score > out[i]:
                    out[i] = score
                    max_clones[i] = set([node[2]])
                elif score == out[i]:
                    max_clones[i].add(node[2])
                else:
                    break
    return out, max_clones

class ScoreTimeline(Analyzer):
    def __init__(self, project_name, test, num_timepoints=100):
        super().__init__(project_name, test)
        self.num_timepoints = num_timepoints

    def analyze(self):
        return get_max_by_time(self.cached_trees[-1], self.num_timepoints)