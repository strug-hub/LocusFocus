/**
 * Populates the #set-based-test-table with results from the first stage significance test.
 */
function buildFirstStageTable(sessionData) {
    let testResultData = [];
    let descriptionHeader = "";
    if (sessionData.hasOwnProperty('regions')) {
        // set-based only
        const numTests = sessionData["first_stages"].length;
        const numRegions = sessionData["regions"].length;
        const multipleTestsRequested = sessionData["multiple_tests"]
        if (multipleTestsRequested) {
            // multiple tests
            descriptionHeader = "Test Region";
            testResultData = sessionData["regions"].map((regiontext, i) => [
                regiontext,
                sessionData["first_stages"][i] ? "Yes" : "No",
                sessionData["first_stage_Pvalues"][i],
            ])
        } else {
            // should be 1 test
            if (numTests > 1) console.warn(`'${numTests}' test detected despite not requesting multiple tests.`);
            descriptionHeader = "Dataset description";
            // list of lists
            testResultData = [[
                sessionData["dataset_title"],
                sessionData["first_stages"][0] ? "Yes" : "No",
                sessionData["first_stage_Pvalues"][0],
            ]];
        }        
    } else {
        let titleKey = "secondary_dataset_titles";
        testResultData = sessionData[titleKey].map((title, i) => [
            title,
            sessionData["first_stages"][i] ? "Yes" : "No",
            sessionData["first_stage_Pvalues"][i],
        ]);
    }

    const columns = [
        { title: descriptionHeader },
        { title: "Set-based test passed?" },
        { title: "Set-based test P-value" },
    ];

    $(document).ready(() => {
        $("#set-based-test-table").DataTable({
            dom: "Bfrtipl",
            data: testResultData,
            columns: columns,
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
