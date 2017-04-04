def quantile_plot(x, y, quantiles = 5000, ax = None, s = 5):
    assert x.shape == y.shape
    data_pairs = [(X, Y) for X, Y in zip(x, y)]
    data_pairs = sorted(data_pairs, key = lambda x: x[0])
    avg_activity, avg_score = [], []
    for i in range(0, len(data_pairs), len(data_pairs) / quantiles):
        index = range(i, min(i + (len(data_pairs) / quantiles), len(data_pairs)))
        activities = [data_pairs[j][0] for j in index]
        scores  = [data_pairs[j][1] for j in index]
        avg_activity.append(sum(activities) / float(len(activities)))
        avg_score.append(sum(scores) / float(len(scores)))
    if ax:
        ax.scatter(avg_activity, avg_score, c = 'r', s=s)
    else:
        plt.scatter(avg_activity, avg_score, c = 'r', s=s)
        plt.show()