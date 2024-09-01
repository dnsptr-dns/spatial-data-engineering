CREATE
OR REPLACE VIEW staging.linked_data_view AS
SELECT
    a.ID,
    b."TEMA" TEMA,
    a."LUSE" LUSE,
    a."KETERANGAN" KETERANGAN,
    b."JENIS" JENIS,
    b."SUMBER" SUMBER,
    a.geometry
FROM
    staging.tb_lu_dataset a
    JOIN staging.tb_lu_csv_dataset b ON a."TEMA" = b."TEMA";