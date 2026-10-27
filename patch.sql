-- Sep 25, 2026
-- Introduce two new file types for qtp-genome

INSERT INTO qiita.filepath_type (filepath_type) VALUES ('assembly')
ON CONFLICT (filepath_type) DO NOTHING;

INSERT INTO qiita.filepath_type (filepath_type) VALUES ('annotation')
ON CONFLICT (filepath_type) DO NOTHING;
