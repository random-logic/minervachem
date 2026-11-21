from sklearn.pipeline import Pipeline
from sklearn.linear_model import Ridge
from minervachem.fingerprinters import GraphletFingerprinter
from minervachem.transformers import FingerprintFeaturizer
from sklearn.model_selection import train_test_split

def create_pipeline(df, n_jobs=None):
    """creates pipeline object of the MinervaChem ridge model.

    Args:
        df (pandas dataframe): a dataframe of qm9 molecules that have been processed with the /demos/preprocess_qm9.py

    Returns:
        pipeline: model pipeline object
    """

    if n_jobs is None:
        import os
        n_jobs = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))

    train, test = train_test_split(
        df,
        train_size=0.8,
        random_state=42,
    )

    y_train, _ = [sub_df['E_at'] for sub_df in [train, test]]


    pipeline = Pipeline(
        [
            ("featurizer", FingerprintFeaturizer(fingerprinter=GraphletFingerprinter(max_len=3, useHs=True),
                                                verbose=0,
                                                n_jobs=n_jobs,
                                                chunk_size='auto',)
                                                ),
            ("ridge", Ridge(fit_intercept=False,
                            alpha=1e-5,
                            solver='sparse_cg',
                            tol=1e-5)
                            )
        ]
    )

    pipeline.fit(train['mol'], y_train)

    return pipeline


if __name__ == '__main__':
    import argparse
    from pathlib import Path
    import os
    import joblib
    import pandas as pd
    from rdkit import Chem

    parser = argparse.ArgumentParser(description='Train and save the molecular QM9 ridge pipeline.')
    parser.add_argument('--data', type=Path, default=Path(__file__).resolve().parents[1] / 'qm9_processed.csv')
    parser.add_argument('--output', type=Path, default=Path('qm9_pipeline.joblib'))
    parser.add_argument('--limit', type=int, help='Use the first N molecules for a quick run')
    parser.add_argument('--n-jobs', type=int, default=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    args = parser.parse_args()
    frame = pd.read_csv(args.data)
    if args.limit is not None:
        frame = frame.head(args.limit)
    frame['mol'] = frame['smiles'].map(lambda smiles: Chem.AddHs(Chem.MolFromSmiles(smiles)))
    pipeline = create_pipeline(frame, n_jobs=args.n_jobs)
    joblib.dump(pipeline, args.output)
    print(f'Saved molecular pipeline to {args.output}')
