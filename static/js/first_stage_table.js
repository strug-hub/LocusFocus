/**
 * Populates the #set-based-test-table with results from the first stage significance test.
 */
function buildFirstStageTable(sessionData) {
    let titleKey = "secondary_dataset_titles";
    if (sessionData.hasOwnProperty('dataset_titles')) {
        // for set-based only
        titleKey = 'dataset_titles';
    }

    let testResultData = sessionData[titleKey].map((title, i) => [
        title,
        sessionData["first_stages"][i] ? "Yes" : "No",
        sessionData["first_stage_Pvalues"][i],
    ]);

    $(document).ready(() => {
        $("#set-based-test-table").DataTable({
            dom: "Bfrtipl",
            data: testResultData,
            columns: [
            { title: "Dataset description" },
            { title: "First-stage test passed?" },
            { title: "First-stage test P-value" },
            ],
            buttons: [
            "copy",
            {
                extend: "csv",
                filename: "Secondary_datasets_SS_pvalues",
            },
            {
                extend: "excel",
                filename: "Secondary_datasets_SS_pvalues",
                messageTop: "Simple Sum P-values of Secondary Datasets",
            },
            ],
        });
    });
}

function buildFirstStageGTEXTable(sessionData) {
    // TODO: implement when GTEx is ready for set-based route
    return null;
}