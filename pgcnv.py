import networkx as nx
import torch
import torch.nn as nn
import torch.nn.functional as F
import time
from datetime import datetime
from scipy.sparse import coo_matrix
import pandas as pd
import pickle
import os
import argparse
from utils import *
torch.manual_seed(42)
np.random.seed(42)

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the '.config' file ", required=True)
parser.add_argument('-k', type=str, help="Path to the '.gpickle' file ", required=True)
parser.add_argument('-o', type=str, help="Path to the output file, example: data/output", required=True)

args = parser.parse_args()


class GraphConvolution(nn.Module):
    def __init__(self, in_features, out_features):
        super(GraphConvolution, self).__init__()
        self.linear = nn.Linear(in_features, out_features)

    def forward(self, x, adj):
        x = torch.sparse.mm(adj, x)
        x = self.linear(x)
        return x


class GCN(nn.Module):
    def __init__(self, input_dim, hidden_dim, output_dim):
        super(GCN, self).__init__()
        self.gc1 = GraphConvolution(input_dim, hidden_dim)
        self.gc2 = GraphConvolution(hidden_dim, hidden_dim)
        self.fc = nn.Linear(hidden_dim, output_dim)

    def forward(self, x, adj):
        x = F.relu(self.gc1(x, adj))
        x = F.relu(self.gc2(x, adj))
        x = self.fc(x)
        return x


def sparse_tensor_from_coo(coo):
    values = coo.data
    indices = np.vstack((coo.row, coo.col))
    i = torch.tensor(indices, dtype=torch.long)
    v = torch.tensor(values, dtype=torch.float32)
    shape = coo.shape
    return torch.sparse_coo_tensor(i, v, torch.Size(shape))


def main():
    current_datetime = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    output_path = args.o
    os.makedirs(output_path, exist_ok=True)
    output_file = os.path.join(output_path, f'pgcnv_res_{current_datetime}.tsv')

    with open(args.k, 'rb') as f:
        G = pickle.load(f)

    node_features = []
    labels = []
    nodes = list(G.nodes(data=True))

    for node, data in nodes:
        try:
            feature = [data['start'], data['length'], data['logr']]
            label = data['label']
            node_features.append(feature)
            labels.append(label)
        except KeyError as e:
            print(f"Missing key {e} in node {node}")

    node_features = np.array(node_features, dtype=np.float32)
    label_mapping = {'0': 0, '1': 1, 'unknown': -1}
    labels = np.array([label_mapping[str(label)] for label in labels], dtype=np.float32)

    X = torch.tensor(node_features, dtype=torch.float32)
    y = torch.tensor(labels, dtype=torch.float32)

    adj = nx.adjacency_matrix(G)

    train_idx = np.where(labels != -1)[0]
    test_idx = np.where(labels == -1)[0]

    X_train = X[train_idx]
    y_train = y[train_idx]

    X_test = X[test_idx]

    adj_train = adj[train_idx, :][:, train_idx]
    adj_test = adj[test_idx, :][:, test_idx]

    adj_train = sparse_tensor_from_coo(coo_matrix(adj_train))
    adj_test = sparse_tensor_from_coo(coo_matrix(adj_test))

    input_dim = X.shape[1]
    hidden_dim = 64
    output_dim = 1

    model = GCN(input_dim, hidden_dim, output_dim)

    optimizer = torch.optim.Adam(model.parameters(), lr=float(0.001))
    loss_fn = nn.BCEWithLogitsLoss()  # Binary Cross Entropy Loss

    num_epochs = 200
    for epoch in range(num_epochs):
        model.train()
        optimizer.zero_grad()
        output = model(X_train, adj_train)
        loss = loss_fn(output.squeeze(), y_train)
        loss.backward()
        optimizer.step()
        if (epoch + 1) % 10 == 0:
            print(f"======================Epoch [{epoch + 1}/{num_epochs}]======================")

    model.eval()
    with torch.no_grad():
        output = model(X_test, adj_test)
        predictions = torch.sigmoid(output).squeeze().numpy()

    torch.save(model.state_dict(), 'model.pth')
    print("Model saved as 'model.pth'")

    predicted_labels = np.round(predictions).astype(int)
    test_nodes = np.array(nodes)[test_idx]

    results = []
    for i, (node, data) in enumerate(test_nodes):
        result = {
            'SampleID': node,
            'Chromosome': data['chr_name'],
            'Start': data['start'],
            'End': data['start'] + data['length'],
            'LogR_Ratio': data['logr'],
            'Predicted_Label': predicted_labels[i]
        }
        results.append(result)

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results have saved to the '{output_file}' file")


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")
