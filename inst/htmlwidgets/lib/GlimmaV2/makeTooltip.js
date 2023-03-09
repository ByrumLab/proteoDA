function makeVegaTooltip(columns)
{
    // generate tooltip object for embedding in spec
    let tooltipString = "{";
    columns.forEach( x => {if (x !== "internal_id_for_brushing") {tooltipString += `'${x}':datum['${x}'],`}});
    tooltipString += "}";
    var tooltip = { "signal" : tooltipString };
    return tooltip;
}
