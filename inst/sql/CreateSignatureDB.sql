-- !preview conn=DBI::dbConnect(RSQLite::SQLite())

-- Decompositions
CREATE TABLE decompositions (
    sampleId TEXT,
    class TEXT,
    channel TEXT,
    fraction FLOAT,
    count FLOAT,
    unique (sampleId, channel, class)
);

CREATE INDEX idx_sample_class_channel ON decompositions (sampleId, class, channel);
CREATE INDEX idx_class ON decompositions (class);
CREATE INDEX idx_fraction_count ON decompositions (fraction, count);
CREATE INDEX idx_count_fraction ON decompositions (count, fraction);


-- Cosmic Exposures
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

-- error & cosine
-- "Method","SampleID","Type","IsOptimal","Errors", "Cosine"
CREATE TABLE cosmicErrorAndCosine (
    sampleId TEXT,
    class TEXT,
    type FLOAT,
    optimal BOOLEAN,
    method TEXT CHECK (method IN ('SA', 'QP')),
    error FLOAT,
    cosine FLOAT,
    unique (sampleId, type, method, class)
);

CREATE INDEX idx_sample_optimal ON cosmicErrorAndCosine (sampleId, optimal);
CREATE INDEX idx_optimal_sample ON cosmicErrorAndCosine (optimal, sampleId);
CREATE INDEX idx_optimal_error ON cosmicErrorAndCosine (optimal, error);
CREATE INDEX idx_optimal_cosine ON cosmicErrorAndCosine (optimal, cosine);




-- p_val
CREATE TABLE cosmicPvalues (
    sampleId TEXT,
    class TEXT,
    signature TEXT,
    method TEXT CHECK (method IN ('SA', 'QP')),
    threshold FLOAT,
    p FLOAT,
    unique (sampleId, signature, method, class)
);


CREATE INDEX idx_p_sample_method ON cosmicPvalues (p, sampleId, method);
CREATE INDEX idx_p_signature_method ON cosmicPvalues (p, signature, method);
CREATE INDEX idx_signature_p_method ON cosmicPvalues (signature, p, method);
CREATE INDEX idx_signature_method ON cosmicPvalues (signature, p, method);
CREATE INDEX idx_sample_signature_method ON cosmicPvalues (sampleId, signature, method);

-- Pairwise Cosine Similarity
CREATE TABLE pairwiseSimilarity (
    sample1 TEXT,
    sample2 TEXT,
    class TEXT,
    cosine_similarity FLOAT,
    unique (sample1, sample2, class),
    PRIMARY KEY (sample1, sample2, class)
);

CREATE INDEX idx_pairwise_sample1cosine ON pairwiseSimilarity (sample1, cosine_similarity);
CREATE INDEX idx_pairwise_sample2cosine ON pairwiseSimilarity (sample2, cosine_similarity);
CREATE INDEX idx_pairwise_sample1and2 ON pairwiseSimilarity (sample1, sample2);
CREATE INDEX idx_pairwise_sample2and1 ON pairwiseSimilarity (sample2, sample1);
CREATE INDEX idx_pairwise_cosine_class ON pairwiseSimilarity (cosine_similarity, class);


-- Add a sample-level view describing total mutation count for each decomposition
CREATE VIEW sample AS
  SELECT
      sampleId,
      class,
      SUM(count) AS total_count
  FROM
      decompositions
  GROUP BY
      sampleId,
      class;

