HTMLWidgets.widget({

  name: 'glimmaXY',

  type: 'output',

  factory: function(el, width, height)
  {

    var plotContainer = document.createElement("div");
    var controlContainer = document.createElement("div");
    plotContainer.setAttribute("class", "plotContainer");
    controlContainer.setAttribute("class", "controlContainer");

    var widget = document.getElementById(el.id);
    widget.appendChild(plotContainer);
    widget.appendChild(controlContainer);

    return {

      renderValue: function(x)
      {

        //console.log(x);
        var handler = new vegaTooltip.Handler();

        // create container elements
        var xyContainer = document.createElement("div");
        xyContainer.setAttribute("class", "xyContainerSingle");
        plotContainer.appendChild(xyContainer);

        var xyTable = HTMLWidgets.dataframeToD3(x.data.table)
        var xySpec = createXYSpec(x.data, xyTable, width, height);
        var xyView = new vega.View(vega.parse(xySpec), {
          renderer: 'svg',
          container: xyContainer,
          bind: controlContainer,
          hover: true
        });
        xyView.tooltip(handler.call);
        xyView.runAsync();

        var countsMatrix = null;
        var expressionView = null;
        var expressionContainer = null;
        if (x.data.counts != -1)
        {
          expressionContainer = document.createElement("div");
          expressionContainer.setAttribute("class", "expressionContainer");
          plotContainer.appendChild(expressionContainer);
          xyContainer.setAttribute("class", "xyContainer");
          countsMatrix = HTMLWidgets.dataframeToD3(x.data.counts);
          var expressionSpec = createExpressionSpec(width, height, x.data.expCols, x.data.sampleColours, x.data.samples);
          var expressionView = new vega.View(vega.parse(expressionSpec), {
            renderer: 'svg',
            container: expressionContainer,
            hover: true
          });
          expressionView.tooltip(handler.call);
          expressionView.runAsync();
        }

        var data =
        {
          xyView: xyView,
          expressionView: expressionView,
          xyTable: xyTable,
          countsMatrix: countsMatrix,
          controlContainer: controlContainer,
          height: height,
          cols: x.data.cols,
          groups: x.data.groups,
          levels: x.data.levels,
          expressionContainer: expressionContainer,
          width: width
        };

        setupXYInteraction(data);
        addSavePlotButton(controlContainer, xyView, expressionView, "Save Plot");
        if (expressionView) {
          addAxisMessage(data);
        }
      },

      resize: function(width, height)
      {}

    };
  }
});

class State {

  /**
   * Returns state machine object retaining the current set of selected genes and managing
   * whether the app is in graph selection mode or table selection mode
   * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
   * @return {State} state machine object
   */
  constructor(data) {
    this.data = data;
    this.graphMode = false;
    this._selected = [];
  }

  /**
   * Returns current selection of genes
   * @return {Array} Array of currently selected genes
   */
  get selected() {
    return this._selected;
  }

  /**
   * Sets a new array of selected genes and re-renders elements accordingly
   * @param  {Array} selected Array of genes which are currently selected
   */
  set selected(selected) {
    this._selected = selected;
    let htmlString = selected.map(x => `<span>${x.gene}</span>`).join("");
    $(this.data.controlContainer.getElementsByClassName("geneDisplay")[0])
      .html(htmlString);
    /* update save btn */
    $(this.data.controlContainer.getElementsByClassName("saveSelectButton")[0])
      .html(`Save (${selected.length})`);
    /* update clear btn */
    $(this.data.controlContainer.getElementsByClassName("clearSubset")[0])
      .html(`Clear (${selected.length})`);
  }

  /**
   * Adds a gene to the selection if it's not already selected, or remove it otherwise
   * @param  {Gene} gene Gene data object which has been clicked on
   */
  toggleGene(gene) {
    let loc = containsGene(this.selected, gene);
    this.selected = loc >= 0 ? remove(this.selected, loc) : this.selected.concat(gene);
    this._expressionUpdateHandler(loc < 0, gene);
  }

  /**
   * Manages updates to the expression plot based on the most recently selected gene
   * @param {Boolean} selectionOccurred True if a gene was selected, false if it was de-selected
   * @param  {Gene} gene Gene data object which has been clicked on
   */
  _expressionUpdateHandler(selectionOccurred, gene) {
    if (!this.data.expressionView) return;
    if (selectionOccurred) {
      let countsRow = this.data.countsMatrix[gene.index];
      updateExpressionPlot(countsRow, this.data, gene.gene);
    }
    else if (this.selected.length > 0) {
      let last = this.selected[this.selected.length-1];
      let countsRow = this.data.countsMatrix[last.index];
      updateExpressionPlot(countsRow, this.data, last.gene);
    }
    else {
      clearExpressionPlot(this.data);
    }
  }

}

/**
 * Generates datatable DOM object, state machine and assigns event listeners
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function setupXYInteraction(data)
{

  var state = new State(data);
  var datatableEl = document.createElement("TABLE");
  datatableEl.setAttribute("class", "dataTable");
  data.controlContainer.appendChild(datatableEl);

  $(document).ready(function()
  {
    var datatable = $(datatableEl).DataTable(
      {
        data: data.xyTable,
        columns: data.cols.map(el => ({"data": el, "title": el})),
        columnDefs: [
          {
            targets: 0,
            visible: false,
          },
        ],
        rowId: "gene",
        dom: '<"geneDisplay fade-in">Bfrtip',
        buttons: {
          dom: {
            buttonContainer: {
              tag: 'div',
              className: 'buttonContainer'
            }
          },
          buttons: [
                    {
                      text: 'Clear (0)',
                      action: () => clearTableListener(datatable, state, data),
                      attr: {class: 'save-button clearSubset'}
                    },
                    {
                      text: 'Save Data',
                      action: () => showDataDropdown(),
                      attr: {class: 'save-button saveSubset'}
                    }
                  ]
                },
        scrollY: (data.height*0.4).toString() + "px",
        scrollX: false,
        orderClasses: false,
        stripeClasses: ['stripe1','stripe2']
      });

    var col_number = datatable.columns().count();
    // Loop of column indices, get column names
    for (var i = 0; i < col_number; i++) {
      var title = $(datatable.column(i).header()).text();
      if (title === "negLog10adjP") {
        $(datatable.column(i).visible(false));
      }
      if (title === "negLog10rawP") {
        $(datatable.column(i).visible(false));
      }
    }

    datatable.on('click', 'tr', function() { tableClickListener(datatable, state, data, $(this)) } );
    data.xyView.addSignalListener('click', function(name, value) { XYSignalListener(datatable, state, value[0], data) } );

    $(document.getElementsByClassName("saveSubset")[0]).html(`Save Data`);
    addSaveDataElement(state, data, `Save All`, `Save (0)`);
  });
}

/**
 * Shows Save Data options
 */
function showDataDropdown() {
  let dataDropdown = document.getElementsByClassName("dataDropdown")[0];
  dropdownOnClick(dataDropdown);
}

/**
 * Responds to a click on the Clear datatable button
 * @param  {Datatable} datatable datatable object
 * @param  {State} state state machine object returned by getStateMachine()
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function clearTableListener(datatable, state, data)
{
  state.graphMode = false;
  state.selected = [];
  datatable.rows('.selected').nodes().to$().removeClass('selected');
  datatable.search('').columns().search('').draw();
  data.xyView.data("selected_points", state.selected);
  data.xyView.runAsync();
  clearExpressionPlot(data);
  //console.log(state);
}

/**
 * Listens and responds to click events on the datatable
 * @param  {Datatable} datatable datatable object
 * @param  {State} state state machine object returned by getStateMachine()
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 * @param  {Row} row row object in the table clicked on by the user
 */
function tableClickListener(datatable, state, data, row)
{
  if (state.graphMode) return;
  row.toggleClass('selected');
  let datum = datatable.row(row).data();
  state.toggleGene(datum);
  data.xyView.data("selected_points", state.selected);
  data.xyView.runAsync();
}

/**
 * Listens and responds to click events on the XY plot
 * @param  {Datatable} datatable datatable object
 * @param  {State} state state machine object returned by getStateMachine()
 * @param  {Datum} datum point on the graph clicked on by the user
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function XYSignalListener(datatable, state, datum, data)
{
  if (datum == null) return;
  if (!state.graphMode)
  {
    state.graphMode = true;
    datatable.rows('.selected').nodes().to$().removeClass('selected');
    state.selected = [];
  }

  state.toggleGene(datum);

  // edge case: deselecting last point
  if (state.selected.length == 0)
    state.graphMode = false;

  data.xyView.data("selected_points", state.selected);
  data.xyView.runAsync();

  datatable.search('').columns().search('').draw();
  var regex_search = state.selected.map(x => '^' + x.gene + '$').join('|');
  datatable.columns(0).search(regex_search, regex=true, smart=false).draw();
}

/**
 * Resets expression plot to a blank slate
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function clearExpressionPlot(data)
{

  if (!data.expressionView)
    return;

  data.expressionView.data("table", []);
  data.expressionView.signal("title_signal", "");
  data.expressionView.signal("max_count", 0);
  data.expressionView.runAsync();
  updateAxisMessage(data);
}

/**
 * Updates expression plot for the given gene and sample counts
 * @param  {CountsRow} countsRow Data object containing sample counts for a given gene
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 * @param  {String} geneName name of gene being displayed
 */
function updateExpressionPlot(countsRow, data, geneName)
{
  let groups = data.groups.group;
  let numUniqueGroups = [...new Set(groups)].length;
  let samples = data.groups.sample;
  let levels = data.levels;
  let result = [];
  for (col in countsRow)
  {
    if (!samples.includes(col)) continue;
    let curr = {};
    let group = groups[samples.indexOf(col)];
    curr["group"] = group;
    curr["sample"] = col;
    curr["normalized intensity"] = countsRow[col];
    // Get an offset for jittering
    // First, randomly sample the sign
    Math.random() < 0.5 ? sign = -1 : sign = 1;
    // Then, calculate a jitter width scaling.
    // Vega is using point scale with a padding of 1,
    // So space between points is the total pixel width
    // divided by number of groups + 1
    // So, the maximum we can jitter in wither direction is half of that.
    // Will say we can jitter within 60% of that max right now.
    // But might want to adjust that by number of groups: e.g., jitter even less
    // when there are more groups.
    jit_width = (data.width*0.4)/(numUniqueGroups + 1)/2*0.60;
    // Then, get the offset. Pretty simple: get a random number between 0 and 1, scale
    // by width, and then multiple by the randomly generated sign
    curr["offset"] = Math.random()*jit_width*sign;
    if (!Object.values(curr).includes(-9)) { // if -9, which we're using as missing, don't push result
       result.push(curr);
    }
  }
  if (levels != null) {
    result.sort((a, b) => levels.indexOf(a.group) - levels.indexOf(b.group));
  }
  // Uncomment next line to log the results to the console for debugging
  //console.log(JSON.parse(JSON.stringify(result)));
  data.expressionView.data("table", result);
  data.expressionView.signal("title_signal", "Gene " + geneName.toString());
  let max_value = Math.max(...result.map(x => x["normalized intensity"]));
  let min_value = Math.min(...result.map(x => x["normalized intensity"]));
  data.expressionView.signal("max_count", Math.round(max_value*100)/100 );
  data.expressionView.signal("min_count", Math.round(min_value*100)/100 );
  data.expressionView.runAsync();
  updateAxisMessage(data);
}

/**
 * Adds y-axis scaling message DOM objects to the expression plot
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function addAxisMessage(data)
{
  var bindings = data.expressionContainer.getElementsByClassName("vega-bindings")[0];
  var alertBox = document.createElement("div");
  alertBox.setAttribute("class", "alertBox invisible");
  data.expressionView.addSignalListener('ylim_max',
    function(name, value) { updateAxisMessage(data) });
  data.expressionView.addSignalListener('ylim_min',
    function(name, value) { updateAxisMessage(data) });
  bindings.appendChild(alertBox);
}

/**
 * Updaes the y-axis scaling for the expression plot
 * @param  {Data} data encapsulated data object containing references to Vega graphs and DOM elements
 */
function updateAxisMessage(data)
{
  var alertBox = data.expressionContainer.getElementsByClassName("alertBox")[0];
  let minCount = data.expressionView.signal("min_count");
  let maxCount = data.expressionView.signal("max_count");
  let userMax = data.expressionView.signal("ylim_max");
  let userMin = data.expressionView.signal("ylim_min");

  max_OK = userMax == null || userMax == "" || Number(userMax) >= maxCount;
  min_OK = userMin == null || userMin == "" || Number(userMin) <= minCount;
  if (max_OK && min_OK)
  {
    alertBox.setAttribute("class", "alertBox invisible");
  }
  else
  {
    if (!max_OK)
    {
      alertBox.innerHTML = `Max normalized intensity is ${maxCount}`;
    }
    else
    {
      alertBox.innerHTML = `Min normalized intensity is ${minCount}`;
    }
    alertBox.setAttribute("class", "alertBox danger");
  }
}

/**
 * Searches an array gene data objects to determine if it contains a given gene.
 * @param  {Array} arr array of gene data objects.
 * @param  {Datum} datum given gene object
 * @return {Integer} -1 if the given gene is not found; index of the gene in arr otherwise.
 */
function containsGene(arr, datum)
{
  let loc = -1;
  let i;
  for (i = 0; i < arr.length; i++)
  {
    if (arr[i]['gene'] === datum['gene'])
    {
      loc = i;
      break;
    }
  }
  return loc;
}

/**
 * Removes an element at the given index from an array and returns the result.
 * @param  {Array} arr array of elements.
 * @param  {Integer} i index i of element to be removed from arr.
 * @return {Array} modified array with element at index i removed.
 */
function remove(arr, i)
{
  let new_arr = arr.slice(0, i).concat(arr.slice(i+1))
  return new_arr;
}
