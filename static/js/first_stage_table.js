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
        columns.push({ 
            title: "Number of SNPs used", 
            // className: "snp_modal_button", 
            createdCell: (cell, cellData, rowData, rowIndex, colIndex) => {
                $(cell)
                .addClass("snp_modal_button")
                .attr("data-toggle", "modal")
                .attr("data-target", `#snp_modal_${rowIndex}`);
            }
        });
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

        if (multipleTestsRequested) {
            sessionData["regions"].forEach((region, i) => {
                _createSNPModal(i, region, positions[i]);
            })
        }
    });
}

function _createSNPModal(i, region, positions) {
    // create a modal with a datatable inside it
    // there's so many parts that it's just easier to write it as a string...
    const modal = $.parseHTML(
`
<div class="modal fade" id="snp_modal_${i}" tabindex="-1" role="dialog" aria-hidden="true">
  <div class="modal-dialog modal-lg" role="document">
    <div class="modal-content">
      <div class="modal-header">
        <h4 class="modal-title" id="snp_modal_${i}_title">SNPs used in region ${region}</h4>
        <button type="button" class="close" data-dismiss="modal" aria-label="Close">
          <span aria-hidden="true">×</span>
        </button>
      </div>
      <div class="modal-body" id="snp_modal_${i}_body">
        <table id="snp_modal_${i}_table"
        class="table table-striped table-bordered table-condensed table-sm sortable"
        cellspacing="0"
        width="90%"
        >
        </table>
      </div>
    </div>
  </div>
</div>
`);
    $('body').append(modal);
    const modalData = positions;
    $(`#snp_modal_${i}_table`).DataTable({
        dom: "Bfrtipl",
        data: modalData,
        columns: [
            { title: "Chromosome" },
            { title: "Basepair Position" }
        ],
        buttons: [
        "copy",
        {
            extend: "csv",
            filename: "SNPs_used_in_set_based_test",
        },
        {
            extend: "excel",
            filename: "SNPs_used_in_set_based_test",
            messageTop: "SNPs used in set based test of Dataset region(s)",
        },
        ],
    });
    // TODO: Insert the modal somewhere in the page, add DataTable inside with all the positions
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
