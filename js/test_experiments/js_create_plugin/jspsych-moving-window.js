/*
 * Example plugin template
 */

jsPsych.plugins["moving-window"] = (function() {

  var plugin = {};

  // parameters used in the plugin
  plugin.info = {
    name: "moving-window",
    parameters: {
      words: {
        type: jsPsych.plugins.parameterType.STRING, // BOOL, STRING, INT, FLOAT, FUNCTION, KEYCODE, SELECT, HTML_STRING, IMAGE, AUDIO, VIDEO, OBJECT, COMPLEX
        default: undefined
      },
      key: {
        type: jsPsych.plugins.parameterType.KEYCODE,
        default: 32 // use spacebar to advance
      }
    }
  }

  // function trial that runs when you call in the plugin
  plugin.trial = function(display_element, trial) {

    // data saving
    var trial_data = {words: trial.words};
    var rt = [];
    var current_word = 0;

    function create_moving_window(words, position){
      //convert the words string into an array with single words in a sentence
      var word_list = words.split(' ');
      var stimulus = word_list.map(function(word, index) {
        if(index==position) {
          return word;
        } else {
          return "-".repeat(word.length);

        }
      }).join(' '); // join all the elements of the array back into a string
      return stimulus
    }

    function show_stimulus(){
      display_element.innerHTML = "<p>" + create_moving_window(trial.words, current_word) + "</p>";

      jsPsych.pluginAPI.getKeyboardResponse({
        callback_function: after_response, //function that is being executed when subject presses a key
        valid_responses: [trial.key], // what keys are allowed to be pressed?
        rt_method: 'performance',
        persist: false,
        allow_held_key: false
      });

    }

    function after_response(response_info){
      rt.push(response_info.rt);
      end_trial();
    }

    // function 
    function end_trial(){

      trial_data.rt = JSON.stringify(rt);
      // end trial
      jsPsych.finishTrial(trial_data);

    }

    show_stimulus(); // to start the experiment
    
    
  };

  return plugin;
})();
