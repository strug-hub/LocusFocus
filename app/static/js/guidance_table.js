var tableselect = d3.select('#SSguidance-table');

function buildSSguidanceTable(genes, tissues, SSP, SSP2, smr) {
    var thead = d3.select('#SSguidance-table').select("thead");
    var tbody = d3.select("#SSguidance-table").select("tbody");

    // Clear table:
    if ($.fn.dataTable.isDataTable('#SSguidance-table')) {
        var mytable = $('#SSguidance-table').DataTable();
        mytable.destroy();
    }

    thead.text("");
    tbody.text("");

    // Clear table:
    if ($.fn.dataTable.isDataTable('#SSguidance-table')) {
        var mytable = $('#SSguidance-table').DataTable();
        mytable.destroy();
    }

    // Column headers
    header = thead.append('tr');
    header
        .append('th')
        .attr('class', 'th-sm')
        .text('Field');
    header
        .append('th')
        .attr('class', 'th-sm')
        .text('Value');

    var allSSPs = SSP.flat().concat(smr?.ssp_values || []).concat(SSP2);
    var numGTEx = tissues.length * genes.length;
    var numSMR = smr?.ssp_values?.length || 0;
    var numSecondary = SSP2.length;
    var numNoeQTL = allSSPs.reduce((acc, curr) => (acc + (curr == -1)), 0);
    var numFirstStage = allSSPs.reduce((acc, curr) => (acc + (curr == -2)), 0);
    var numTested = allSSPs.reduce((acc, curr) => (acc + (curr > 0)), 0);
    var numFailed = allSSPs.reduce((acc, curr) => (acc + (curr == -3 || curr == 0)), 0);

    var suggested_SSP = numTested > 0 ? -Math.log10(0.05 / numTested) : "N/A (No datasets were tested successfully)";

    // Table body:
    var row = tbody.append('tr');
    row.append('td').text('Total number of all secondary datasets');
    row.append('td').text(allSSPs.length);
    var row = tbody.append('tr');
    row.append('td').text('Total number of GTEx datasets');
    row.append('td').text(numGTEx);
    var row = tbody.append('tr');
    row.append('td').text('Total number of SMR datasets');
    row.append('td').text(numSMR);
    var row = tbody.append('tr');
    row.append('td').text('Total number of user-uploaded secondary datasets');
    row.append('td').text(numSecondary);
    var row = tbody.append('tr');
    row.append('td').text('Number of datasets with no eQTL data (-1)');
    row.append('td').text(numNoeQTL);
    var row = tbody.append('tr');
    row.append('td').text('Number of datasets not passing first stage (-2)');
    row.append('td').text(numFirstStage);
    var row = tbody.append('tr');
    row.append('td').text('Number of datasets with computation error (-3)');
    row.append('td').text(numFailed);
    var row = tbody.append('tr');
    row.append('td').text('Number of datasets tested for colocalization');
    row.append('td').text(numTested);
    var row = tbody.append('tr');
    row.append('td').text('Suggested Simple Sum colocalization threshold at alpha 0.05 (-log10P)');
    row.append('td').text(suggested_SSP);


    // Add DataTables functionality:
    SSguidanceTable = $(document).ready(function () {
        var thedatatable3 = $('#SSguidance-table').DataTable({
            dom: 'Bfrtipl',
            buttons: [
                'copy',
                {
                    extend: 'csv',
                    filename: 'SSguidanceTable'
                },
                {
                    extend: 'excel',
                    filename: 'SSguidanceTable',
                    messageTop: 'Guidance table SS computation'
                }
            ]
        });
    });
}

