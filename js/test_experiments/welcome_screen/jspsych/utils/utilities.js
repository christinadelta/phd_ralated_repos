// create function for saving participant responses

function downloadCSV(csv, filename) {
    var csvFile;
    var downloadLink; // container for the document

    // retrieve csv file from experiment
    csvFile = new Blob([csv], {type: "text/csv"});

    //download link
    downloadLink = document.createElement("a");

    // retrieve filename:
    downloadLink.download = filename;

    // create a link to the file
    downloadLink.href = window.URL.createObjectURL(csvFile);

    // hide download link
    downloadLink.style.display = 'none';
    document.body.appendChild(downloadLink);
    downloadLink.click();
}