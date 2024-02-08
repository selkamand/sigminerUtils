-- !preview conn=DBI::dbConnect(RSQLite::SQLite())


-- Sample Level Metadata
CREATE TABLE sample (
    sampleId TEXT NOT NULL,
    disease TEXT NOT NULL,
    description TEXT,
    PRIMARY KEY (sampleId),
    unique(sampleId)
);

CREATE INDEX idx_samplemeta_sample ON sample (sampleId);
CREATE INDEX idx_samplemeta_disease ON sample (disease);

-- Decompositions
CREATE TABLE decompositions (
    sampleId TEXT NOT NULL,
    class TEXT NOT NULL,
    channel TEXT NOT NULL,
    fraction FLOAT,
    count FLOAT NOT NULL,
    unique (sampleId, channel, class)
    PRIMARY KEY (sampleId, channel, class)
    FOREIGN KEY(sampleId) REFERENCES sample(sampleId)
);

CREATE INDEX idx_sample_class_channel ON decompositions (sampleId, class, channel);
CREATE INDEX idx_class ON decompositions (class);
CREATE INDEX idx_fraction_count ON decompositions (fraction, count);
CREATE INDEX idx_count_fraction ON decompositions (count, fraction);


-- Cosmic Exposures
CREATE TABLE cosmicExposures (
    sampleId TEXT NOT NULL,
    signature TEXT NOT NULL,
    contribution FLOAT,
    contributionRelative FLOAT,
    type FLOAT,
    optimal BOOLEAN NOT NULL,
    method TEXT CHECK (method IN ('SA', 'QP')) NOT NULL,
    unique (sampleId, signature, type, method)
    PRIMARY KEY (sampleId, type, signature, method)
    FOREIGN KEY(sampleId) REFERENCES sample(sampleId)
);

CREATE INDEX idx_sample_signature_optimal ON cosmicExposures (sampleId, signature, optimal);
CREATE INDEX idx_signature_optimal ON cosmicExposures (signature, optimal);
CREATE INDEX idx_optimal_contributionRelative_signature ON cosmicExposures (optimal, contributionRelative, Signature);

-- error & cosine
-- "Method","SampleID","Type","IsOptimal","Errors", "Cosine"
CREATE TABLE cosmicErrorAndCosine (
    sampleId TEXT NOT NULL,
    class TEXT NOT NULL,
    type FLOAT NOT NULL,
    optimal BOOLEAN NOT NULL,
    method TEXT CHECK (method IN ('SA', 'QP')) NOT NULL,
    error FLOAT,
    cosine FLOAT,
    unique (sampleId, type, method, class)
    FOREIGN KEY(sampleId) REFERENCES sample(sampleId)
);

CREATE INDEX idx_sample_optimal ON cosmicErrorAndCosine (sampleId, optimal);
CREATE INDEX idx_optimal_sample ON cosmicErrorAndCosine (optimal, sampleId);
CREATE INDEX idx_optimal_error ON cosmicErrorAndCosine (optimal, error);
CREATE INDEX idx_optimal_cosine ON cosmicErrorAndCosine (optimal, cosine);




-- p_val
CREATE TABLE cosmicPvalues (
    sampleId TEXT NOT NULL,
    class TEXT NOT NULL,
    signature TEXT NOT NULL,
    method TEXT CHECK (method IN ('SA', 'QP')) NOT NULL,
    threshold FLOAT NOT NULL,
    p FLOAT,
    unique (sampleId, signature, method, class)
    FOREIGN KEY(sampleId) REFERENCES sample(sampleId)
);


CREATE INDEX idx_p_sample_method ON cosmicPvalues (p, sampleId, method);
CREATE INDEX idx_p_signature_method ON cosmicPvalues (p, signature, method);
CREATE INDEX idx_signature_p_method ON cosmicPvalues (signature, p, method);
CREATE INDEX idx_signature_method ON cosmicPvalues (signature, p, method);
CREATE INDEX idx_sample_signature_method ON cosmicPvalues (sampleId, signature, method);

-- Pairwise Cosine Similarity
CREATE TABLE pairwiseSimilarity (
    sample1 TEXT NOT NULL,
    sample2 TEXT NOT NULL,
    class TEXT NOT NULL,
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
CREATE VIEW sample_class_summaries AS
  SELECT
      sampleId,
      class,
      SUM(count) AS total_count
  FROM
      decompositions
  GROUP BY
      sampleId,
      class;

