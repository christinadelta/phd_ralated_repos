// Experimental parameters are defined here

// CONSTANTS
var TODAY = new Date();
var DD = String(TODAY.getDate()).padStart(2, '0');
var MM = String(TODAY.getMonth() + 1).padStart(2, '0');
var YYYY = TODAY.getFullYear();
const DATE = DD + MM + YYYY;

//experimental params
PRACT_REP = 2;
EXP_REP = 2;

// 1. welcome screen
var welcome_block = {
    data: {
        screen_id: "survey-html-form",
    },
    type: "survey-html-form",
    preamble: "<p>Welcome to the experiment!</p>" +
        "Please complete the form",
    html: "<p>Participant ID: <input name='Part_ID' type='text' /></p>",
    on_finish: function (data) {
        console.log(full_design);
        data.responses = JSON.parse(data.responses);
        jsPsych.data.addProperties({Part_ID: data.responses.Part_ID})
    }
};

// init timeline
timeline = [];

// 2. instructions for practice trials
var instructions_block = {
    data: {
        screen_id: "instructions"
    },
    type: "instructions",
    pages: [
        // page 1
        "<p>In this experiment a circle will appear in the middle of the screen.</p>" +
        "<p>If the colour of the circle is <b>blue</b> press <b>left</b> key. </p>" +
        "<p>If the colour of the circle is <b>orange</b> press the <b>right</b> key. </p>" +
        "<div style='float: left';> <img src='jspsych/examples/img/blue.png'/>" +
        "<p><b>Press the left key</b></p></div>" +
        "<div style='float: right';> <img src='jspsych/examples/img/orange.png'/>" +
        "<p><b>Press the right key<b></b></p></div>",

        // page 2
        "We will begin with some practice trials."
    ],
    show_clickable_nav: true
};

// set trials -- put all stimuli into an array and obtain the stimuli from that array
/* var all_trial_stimuli = [
    {stimulus: "jspsych/examples/img/blue.png", data:{screen_id: "trial", stimulus: "blue", correct_response: 37}},
    {stimulus: "jspsych/examples/img/blue.png", data:{screen_id: "trial", stimulus: "blue", correct_response: 37}},
    {stimulus: "jspsych/examples/img/orange.png", data:{screen_id: "trial", stimulus: "orange", correct_response: 39}},
    {stimulus: "jspsych/examples/img/orange.png", data:{screen_id: "trial", stimulus: "orange", correct_response: 39}}
]; */

// set factorial design
factors = {
    stimulus: ["jspsych/examples/img/blue.png", "jspsych/examples/img/orange.png"], // factor 1
    trial_duration: [400, 1000], // factor 2
    // fixation_duration: [400, 1000], // factor 3

};

// full design
var full_design = jsPsych.randomization.factorial(factors, 1);

// loop through trials
var i;
for(i = 0; i < full_design.length; i++) {
    // if blue then data.stimulus: blue, correct responce: 37
    if (full_design[i].stimulus == "jspsych/examples/img/blue.png"){
        full_design[i].data = {screen_id: "trial", stimulus: "blue", correct_response: 37};
    } else { // if blue then data.stimulus: orange, correct responce: 39
        full_design[i].data = {screen_id: "trial", stimulus: "orange", correct_response: 39};
    }
}


// add a fixation point
var fixation = {
    data: {screen_id: "fixation", stimulus: "fixation"},
    type: "html-keyboard-response",
    stimulus: "<div style='font-size:60px'><b>+</b></div>",
    choices: jsPsych.NO_KEYS,
    trial_duration: 500 // one second
};

// 3. set of practice trials
/*
var blue_trial = {
    type: 'image-keyboard-response',
    stimulus: "jspsych/examples/img/blue.png",
    choices: [37,39],
    trial_duration: null
};

var orange_trial = {
    type: 'image-keyboard-response',
    stimulus: "jspsych/examples/img/orange.png",
    choices: [37,39],
    trial_duration: null
}; */

// create the trial variable (includes both blue and orange trials
var trial = {
    data: jsPsych.timelineVariable("data"),
    type: 'image-keyboard-response',
    stimulus: jsPsych.timelineVariable("stimulus"),
    choices: [37,39],
    trial_duration: jsPsych.timelineVariable("trial_duration"),
    on_finish: function(data) {
        //define correct or incorrect (match coorect response with key_press)
        if (data.correct_response == data.key_press) {
            data.accuracy = 1
        } else {
            data.accuracy = 0
        }
    }
};

// create a variable for jittered ITI duration
var ITI_jitt = [500, 750, 1000];

// create feedback screen
var feedback_trial = {
    on_start: function(trial) {
        //get response info from last trial:
        var last_trial = jsPsych.data.get().last(1).values()[0].accuracy
        if (last_trial == 1) { // if accuracy is 1, response was correct -- congratulate participant
            var feedback_text = "Correct!"
        } else {  // if accuracy is 0, response was incorrect -- show incorrect screen
            var feedback_text = "Incorrect!"
        }

        var fback_trial = "<div style='font'size: 90px;'><b>" + feedback_text + "</b></div>";

        trial.data = {screen_id: "feedback", stimulus: feedback_text};
        trial.stimulus = fback_trial;
    },
    data: "",
    type: "html-keyboard-response",
    stimulus: "",
    choices: jsPsych.NO_KEYS,
    trial_duration: 1000,  // 5 seconds
    post_trial_gap: jsPsych.randomization.sampleWithoutReplacement(ITI_jitt, 1) // it will choose one of the values of the ITI_jitt variable
};
// create the order of stimulus presentation:
var procedure = {
    timeline: [fixation, trial, feedback_trial],
    timeline_variables: full_design,
    randomize_order: true, // randomise trials
    repetitions: PRACT_REP // this should be 5 repetition of each colour/condition
};

// experimental instructions ( this is very similar to practice instructions)
var exp_instructions = {
    data: {
        screen_id: "exp_instructions"
    },
    type: "instructions",
    pages: [
        // page 1
        "<p>We have completed our practice trials. Remember,</p>" +
        "<p>If the colour of the circle is <b>blue</b> press <b>left</b> key. </p>" +
        "<p>If the colour of the circle is <b>orange</b> press the <b>right</b> key. </p>" +
        "<div style='float: left';> <img src='jspsych/examples/img/blue.png'/>" +
        "<p><b>Press the left key</b></p></div>" +
        "<div style='float: right';> <img src='jspsych/examples/img/orange.png'/>" +
        "<p><b>Press the right key<b></b></p></div>",

        // page 2
        "We will begin experiment soon."
    ],
    show_clickable_nav: true
};

// create experimental pocedure and trials
var exp_procedure = {
    timeline: [fixation, trial],
    timeline_variables: full_design,
    randomize_order: true, // randomise trials
    repetitions: EXP_REP // this should be 5 repetition of each colour/condition
};

// give debriefing
var debrief = {
    data: {
        screen_id: "debrief"
    },
    type: "instructions",
    pages: [
        // page 1
        "<p>We have completed the experiment. Remember,</p>" +
        "<p>you need to close the screen. You can now contact the experimenter. </p>"
    ],
    show_clickable_nav: true
};

// push to timeline
timeline.push(welcome_block);
timeline.push(instructions_block);
timeline.push(procedure);
timeline.push(exp_instructions);
timeline.push(exp_procedure);
timeline.push(debrief);