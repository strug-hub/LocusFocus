/**
 * Handle checking the status of a job and
 * updating the progress bar on the waiting page.
 */
async function handleJobStatus(jobStatusURL, sessionId) {
  let checks = 0;
  let jobStatus = "PENDING";
  let redirectUrl = "";

  // Wait 2 seconds before checking the status
  await new Promise((resolve) => setTimeout(resolve, 2000));

  // Pending loop
  while (jobStatus == "PENDING") {
    let response = await fetch(jobStatusURL);
    let data = await response.json();
    jobStatus = data.status;
    if (jobStatus == "PENDING") {
      // wait 10 seconds
      await new Promise((resolve) => setTimeout(resolve, 10000));
    }
  }

  // Running loop
  while (jobStatus == "RUNNING") {
    checks += 1;
    let response = await fetch(jobStatusURL);
    let data = await response.json();
    jobStatus = data.status;
    redirectUrl = data.redirect_url ? data.redirect_url : "";
    if (jobStatus == "RUNNING") {
      let percent = 100 - (1 / (1 + (checks + Math.random()) * 0.1)) * 100;
      if (percent > 95) {
        percent = 95;
      }
      document.getElementById("progress-bar").style.width = percent + "%";
      document.getElementById("progress-bar").ariaValueNow = percent;
      // wait 10 seconds
      await new Promise((resolve) => setTimeout(resolve, 10000));
    }
  };
  // Remove spinner
  document.getElementById("loading-spinner").setAttribute("style", "display:none !important");

  if (jobStatus == "FAILURE") {
    document.getElementById("progress-bar").style.width = "100%";
    document.getElementById("progress-bar").ariaValueNow = 100;
    document.getElementById("progress-bar").classList.add("bg-danger");
    document.getElementById("error-section").style.display = "block";
    document.getElementById("error-title").innerHTML = data.error_title;
    document.getElementById("error-message").innerHTML = data.error_message;
  } else if (jobStatus == "SUCCESS") {
    document.getElementById("progress-bar").style.width = "100%";
    document.getElementById("progress-bar").ariaValueNow = 100;
    document.getElementById("progress-bar").classList.add("bg-success");
    document.getElementById("sessionid").innerHTML = document.getElementById("sessionid").innerHTML.replace(sessionId, `<a href="${redirectUrl}">${sessionId}</a>`);
    document.getElementById("job-status-text").innerHTML = `<i>Your results are ready!</i>`;

    document.getElementById("success-section").style.display = "block";
    document.getElementById("success-button").href = redirectUrl;
  }
}
