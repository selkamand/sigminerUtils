-- !preview conn=DBI::dbConnect(RSQLite::SQLite())

CREATE TABLE cosmicExposures (
    sampleId TEXT,
    signature TEXT,
    contribution FLOAT,
    contributionRelative FLOAT,
    type FLOAT,
    optimal BOOLEAN,
    method TEXT CHECK (method IN ('SA', 'QP')),
    unique (sampleId, signature, type, method)
);

CREATE INDEX idx_sample_signature_optimal ON cosmicExposures (sampleId, signature, optimal);
CREATE INDEX idx_signature_optimal ON cosmicExposures (signature, optimal);
CREATE INDEX idx_optimal_contributionRelative_signature ON cosmicExposures (optimal, contributionRelative, Signature);
