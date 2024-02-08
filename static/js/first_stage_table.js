/**
 * Populates the #set-based-test-table with results from the first stage significance test.
 */
function buildFirstStageTable(sessionData) {
    const isSetBasedTest = sessionData.hasOwnProperty('regions');

    if (isSetBasedTest) {
        return _buildSetBasedTestTable(sessionData);
    } else {
        return _buildColocalizationFirstStageTable(sessionData);
    }
}

function _buildSetBasedTestTable(sessionData) {
    // set-based only
    const numTests = sessionData["first_stages"].length;
    const multipleTestsRequested = sessionData["multiple_tests"];
    const positions = sessionData["snps_used_in_test"];
    let testResultData = [];
    if (multipleTestsRequested) {
        // multiple tests
        descriptionHeader = "Test Region";
        testResultData = sessionData["regions"].map((regiontext, i) => [
            regiontext,
            sessionData["first_stage_Pvalues"][i],
            positions[i].length,
        ]);
    } else {
        // should be 1 test
        if (numTests > 1) console.warn(`'${numTests}' test detected despite not requesting multiple tests.`);
        descriptionHeader = "Dataset description";
        // list of lists
        testResultData = [[
            sessionData["dataset_title"],
            sessionData["first_stage_Pvalues"][0],
        ]];
    }

    let columns = [
        { title: descriptionHeader },
        { title: "Set-based test P-value" },
    ];

    if (multipleTestsRequested) {
        columns.push({ title: "Number of SNPs used" });
    }

    $(document).ready(() => {
        $("#set-based-test-table").DataTable({
            dom: "Bfrtipl",
            data: testResultData,
            columns: columns,
            buttons: [
            "copy",
            {
                extend: "csv",
                filename: "Set_based_test_pvalues",
            },
            {
                extend: "excel",
                filename: "Set_based_test_pvalues",
                messageTop: "Set-based test P-values of Dataset region(s)",
            },
            ],
        });

        if (multipleTestsRequested) {}
    });
}

function _buildColocalizationFirstStageTable(sessionData) {
    let titleKey = "secondary_dataset_titles";
    let testResultData = sessionData[titleKey].map((title, i) => [
        title,
        sessionData["first_stage_Pvalues"][i],
    ]);

    let columns = [
        { title: "Secondary dataset" },
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
