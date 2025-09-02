import torch
import joblib
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
def calculate_molecular_features(smiles):
    """Calculate molecular descriptors from SMILES"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    features = [
        Descriptors.MolWt(mol),
        Descriptors.LogP(mol),
        Descriptors.NumHDonors(mol),
        Descriptors.NumHAcceptors(mol),
        # Add more descriptors as needed
    ]
    return np.array(features)
def virtual_screening_pipeline(compound_smiles, model, top_n=50, apply_filters=True, scaler=None):
    """Main screening pipeline"""
    scores = []
    valid_compounds = []

    for smiles in compound_smiles:
        # Calculate molecular features
        features = calculate_molecular_features(smiles)
        if features is None:
            continue

        # Apply drug-like filters (Lipinski's Rule of Five)
        if apply_filters:
            mol = Chem.MolFromSmiles(smiles)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.LogP(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)

            if not (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10):
                continue

        # Scale features if scaler provided
        if scaler:
            features = scaler.transform(features.reshape(1, -1)).flatten()

        # Predict with model
        with torch.no_grad():
            features_tensor = torch.FloatTensor(features).unsqueeze(0)
            score = model(features_tensor).item()

        scores.append(score)
        valid_compounds.append(smiles)

    # Get top hits
    score_indices = np.argsort(scores)[::-1][:top_n]
    top_hits = [(valid_compounds[i], scores[i]) for i in score_indices]

    return top_hits, scores
def analyze_screening_results(top_hits, all_scores):
    """Analyze and visualize screening results"""
    print(f"Screened {len(all_scores)} compounds")
    print(f"Top {len(top_hits)} hits identified")
    print(f"Score range: {min(all_scores):.3f} to {max(all_scores):.3f}")

    print("\nTop 10 compounds:")
    for i, (smiles, score) in enumerate(top_hits[:10]):
        print(f"{i+1}. {smiles} (Score: {score:.3f})") whats s this code saying
