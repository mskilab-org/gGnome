// The configuration parameters
// The Golden Ratio
var phi = 1.618;
var deviation = 2;
var offset = 2;
var slope = 0.25;
var throttleTimer; //used for redrawing upon resize
var totalWidth, totalHeight, width, height, plotsHeight;
var margins = {top: 20, bottom: 50, left: 30, right: 30, gap: 20, bar: 10, legendFontSize: 14, legendHeight: 150, chromosomeContainerHeight: 30, chromobuttonsHeight: 30, chromoAxis: 38, legendGap: 10};
var colorScale = d3.scaleOrdinal(d3.schemeCategory10.concat(d3.schemeCategory20b));
// define the line
var line = d3.line().curve(d3.curveBasis).x(function(d) { return d[0]; }).y(function(d) { return d[1]; });

var metadata = getMetadata(data);
var intervalBins = getIntervalBins(data);
var connectionBins = getConnectionBins(data, intervalBins);
var localConnectionChromosomeBins = getLocalConnectionBins(data, connectionBins);
var looseConnectionChromosomeBins = getLooseConnectionBins(data, connectionBins);
var yMax = d3.max(data.intervals.map(function(d,i) { return d.y; }));
var interChromosomeConnectionBins;
var localInterChromosomeConnectionBins;

// The actual drawing
draw();

d3.select(window).on('resize', throttle);

function draw() {
  // Clear any existing svg
  d3.select('#plot-container svg').remove();

  totalWidth = $('#plot-container').width();
  totalHeight = $(window).height();
  width = totalWidth - margins.left - margins.right;
  height = totalHeight - margins.top - margins.bottom;
  plotsHeight = height - margins.legendHeight - margins.gap;
  
  // The SVG hosting the visualisation
  var svg = d3.select('#plot-container').append('svg').attr('class', 'plot').attr('width', totalWidth).attr('height', totalHeight);

  var panels = data.metadata.slice(0,3).map(function(d,i) { var elem = Object.assign({}, d); elem.column = i; return elem;});

  var panelContainerWidth = (width - (panels.length - 1) * margins.gap) / panels.length;

  drawFilters();

  var yScale = d3.scaleLinear().domain([0, 10, yMax]).range([plotsHeight, 0.4 * plotsHeight, 20]).nice();
  var yAxis = d3.axisLeft(yScale).tickSize(-width).tickValues(d3.range(0, 10).concat(d3.range(10, 10 * Math.round(yMax / 10) + 1, 10)));
  interChromosomeConnectionBins = getInterChromosomeConnectionBins(data, panels, connectionBins);
  localInterChromosomeConnectionBins = getLocalInterChromosomeConnectionBins(data, panels, connectionBins);

  // Controls container
  drawControls();

  // Add the panels
  var panelsContainer = svg.append('g')
    .attr('class', 'panels-container')
    .attr('transform', 'translate(' + [margins.left, height - plotsHeight] + ')');
  
  panelsContainer.append('g')
    .attr('class', 'axis axis--y')
    .attr('transform', 'translate(' + [0, 0] + ')')
    .call(yAxis);

    panelsContainer.select('.axis.axis--y .domain').remove();

  var interChromosomeConnectionsContainer = svg.append('g')
    .attr('class', 'inter-chromosome-connections-container')
    .attr('transform', 'translate(' + [margins.left, height - plotsHeight] + ')');

  updatePanels(panels);

  // Add the legend
  var legendContainer = svg.append('g')
    .attr('class', 'legend-container')
    .attr('transform', 'translate(' + [margins.left, margins.top] + ')');

  updateLegend(panels);

  function drawFilters() {
    var defs = svg.append('defs');

    defs.append('clipPath')
      .attr('id', 'clip')
      .append('rect')
      .attr('width', panelContainerWidth)
      .attr('height', plotsHeight);

    defs.append('clipPath')
      .attr('id', 'clipWidth')
      .append('rect')
      .attr('width', width)
      .attr('height', plotsHeight);

    // create filter and assign provided id
    var filter = defs.append('filter')
      .attr('height', '125%') // adjust this if shadow is clipped
      .attr('id', 'md-shadow');

    // ambient shadow into ambientBlur
    //   may be able to offset and reuse this for cast, unless modified
    filter.append('feGaussianBlur')
      .attr('in', 'SourceAlpha')
      .attr('stdDeviation', deviation)
      .attr('result', 'ambientBlur');

    // cast shadow into castBlur
    filter.append('feGaussianBlur')
      .attr('in', 'SourceAlpha')
      .attr('stdDeviation', deviation)
      .attr('result', 'castBlur');

    // offsetting cast shadow into offsetBlur
    filter.append('feOffset')
      .attr('in', 'castBlur')
      .attr('dx', offset - 1)
      .attr('dy', offset)
      .attr('result', 'offsetBlur');

    // combining ambient and cast shadows
    filter.append('feComposite')
      .attr('in', 'ambientBlur')
      .attr('in2', 'offsetBlur')
      .attr('result', 'compositeShadow');

    // applying alpha and transferring shadow
    filter.append('feComponentTransfer')
      .append('feFuncA')
      .attr('type', 'linear')
      .attr('slope', slope);

    // merging and outputting results
    var feMerge = filter.append('feMerge');
    feMerge.append('feMergeNode')
    feMerge.append('feMergeNode')
      .attr('in', 'SourceGraphic');
  }

  function drawControls() {

    var controlsContainer = svg.append('g')
      .attr('class', 'controls-container')
      .attr('transform', 'translate(' + [margins.left, margins.top] + ')');

    var chromoControls = controlsContainer.selectAll('g.control-container')
      .data(d3.range(panels.length).map(function(d,i) { return {column: i}}), function(d,i) { return d.column })
      .enter()
      .append('g')
      .attr('class', function(d,i) { return 'control-container column-' + i })
      .attr('transform', function(d,i) { return 'translate(' + [i * (panelContainerWidth + margins.gap), 0] + ')'; })

    var chromoButtonContainer = chromoControls
      .append('g')
      .attr('class', 'chromo-buttons-group')
      .selectAll('g.chromo-button-container')
      .data(function(d,i) { return data.metadata.map(function(e,j) {return {chromoObject: e, column: d.column}})}, function(d,i) { return 'column-' + d.column + 'chromo-' + d.chromoObject.chromosome})
      .enter()
      .append('g')
      .attr('class', 'chromo-button-container')
      .attr('transform', function(d,i) { return 'translate(' + [i * (panelContainerWidth -  margins.chromobuttonsHeight)/ (data.metadata.length - 1), ((i + 1) % 2) * margins.chromobuttonsHeight] + ')' });

    chromoButtonContainer.append('circle')
      .attr('class', 'chromo-circle')
      .classed('selected', function(d,i) { return  panels[d.column].chromosome === d.chromoObject.chromosome; })
      .attr('cx', 0.5 * margins.chromobuttonsHeight)
      .attr('cy', 0.5 * margins.chromobuttonsHeight)
      .attr('r', 0.5 * margins.chromobuttonsHeight)
      .attr('fill', function(d,i) { return d.chromoObject.color})
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlight', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlight', false);
      }).on('click', function(d,i) { refreshPanel(d.column, d.chromoObject.chromosome); });

    chromoButtonContainer.append('text')
      .attr('class', 'button-text')
      .attr('transform', 'translate(' + [0.5 * margins.chromobuttonsHeight, 0.5 * margins.chromobuttonsHeight] + ')')
      .attr('text-anchor', 'middle')
      .attr('dy', 0.33 * margins.legendFontSize)
      .text(function(d,i) { return d.chromoObject.chromosome; });
  }

  function updatePanels(newPanels) {

    var panelContainer = panelsContainer.selectAll('g.panel-container').data(newPanels, function(d,i) { return 'column-' + d.column + ' chromo-' + d.chromosome });

    panelContainer.exit().remove();

    var container = panelContainer
      .enter()
      .append('g')
      .attr('class', function(d,i) { return 'panel-container column-' + d.column + ' chromo-' + d.chromosome })
      .attr('transform', function(d,i) { return 'translate(' + [i * (panelContainerWidth + margins.gap), 0] + ')'; })

    container.append('g')
      .attr('class', 'background-container')
      .attr('transform', 'translate(' + [0, 0] + ')')
      .append('rect')
      .attr('class', 'background')
      .attr('width', panelContainerWidth)
      .attr('height', plotsHeight)
      .style('fill', function(d,i) { return metadata[d.chromosome]; });

    container
      .each(function(d,i) {
        d.scale = d3.scaleLinear().domain([metadata[d.chromosome].startPoint, metadata[d.chromosome].endPoint]).range([0, panelContainerWidth]).nice();
        d.panelScale = d3.scaleLinear().domain([metadata[d.chromosome].startPoint, metadata[d.chromosome].endPoint]).range([d.column * (panelContainerWidth + margins.gap), (d.column + 1) * panelContainerWidth + d.column * margins.gap]);
        d.axis = d3.axisBottom(d.scale).ticks(10, 's');
        d.zoom = d3.zoom().scaleExtent([1, Infinity]).translateExtent([[0, 0], [panelContainerWidth, plotsHeight]]).extent([[0, 0], [panelContainerWidth, plotsHeight]]).on('zoom', function() { return zoomed(d)});
        d3.select(this).append('g')
          .attr('class', 'axis axis--x')
          .attr('transform', 'translate(' + [0, plotsHeight] + ')')
          .call(d.axis)
        .selectAll('text')
          .attr('transform', 'rotate(45)')
        .style('text-anchor', 'start');
      });
    
    container.append('rect')
      .attr('class', 'zoom')
      .attr('width', panelContainerWidth)
      .attr('height', plotsHeight)
      .each(function(d,i) {
         d3.select(this).call(d.zoom);
       });

    container.append('g').attr('class', 'shapes-container');

    container.append('g').attr('class', 'local-connections-container');
	
    container.append('g').attr('class', 'local-inter-connections-container');
	
    container.append('g').attr('class', 'loose-connections-container');

  }

  function updateLegend(newPanels) {

    var chromosomeContainer = legendContainer.selectAll('g.chromosome-container').data(newPanels, function(d,i) { return 'column-' + d.column + ' chromo-' + d.chromosome});

    chromosomeContainer.exit().remove();

    var container = chromosomeContainer
      .enter()
      .each(function(d,i) {
        d.scale2 = d3.scaleLinear().domain([metadata[d.chromosome].startPoint, metadata[d.chromosome].endPoint]).range([0, panelContainerWidth]).nice();
        d.axis2 = d3.axisBottom(d.scale2).ticks(10, 's');
      })
      .append('g')
      .attr('class', function(d,i) { return 'chromosome-container column-' + d.column + ' chromosome-' + d.chromosome })
      .attr('transform', function(d,i) { return 'translate(' + [i * (panelContainerWidth + margins.gap), 0] + ')'; })

    container
      .append('g')
      .attr('class', 'axis axis--x')
      .attr('transform', 'translate(' + [0, margins.legendHeight - margins.chromoAxis - 2] + ')')
      .each(function(d,i) { d3.select(this).call(d.axis2).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start'); })

    container
      .append('g')
      .attr('transform', 'translate(' + [0, margins.legendHeight - margins.chromoAxis - margins.chromosomeContainerHeight - margins.legendGap] + ')')
      .append('rect')
      .attr('class', 'chromosome')
      .attr('width', panelContainerWidth)
      .attr('height', margins.chromosomeContainerHeight)
      .style('opacity', function(d,i) { return 0.8; })
      .style('fill', function(d,i) { return d.color; })
      .style('stroke', function(d,i) { return d3.rgb(d.color).darker(1); });
  
    container
      .append('g')
      .attr('transform', 'translate(' + [0, margins.legendHeight - margins.chromoAxis - margins.chromosomeContainerHeight - margins.legendGap] + ')')
      .attr('class', function(d,i) { return 'brush brush-' + d.chromosome; })
      .each(function(d,i) {
        d.brush = d3.brushX().extent([[0, 0], [panelContainerWidth, margins.chromosomeContainerHeight]]).on('brush end', brushed);
        d3.select(this).call(d.brush).call(d.brush.move, d.scale2.range());
      });
  }

  function drawIntervals(panel, scale, dataArray) {

    var shapes = panel.selectAll('rect.shape').data(dataArray, function(d,i) {return d.iid});

    shapes.enter().append('rect')
      .attr('class', 'popovered shape')
      .attr('id', function(d,i) { return 'shape' + d.iid; })
      .style('clip-path','url(#clip)')
      .each(function(d,i) {
        d.startX = scale(d.startPoint);
        d.startY = yScale(d.y);
        d.endX = scale(d.endPoint);
        d.endY = yScale(d.y);
        d.intervalLength = d.endPoint - d.startPoint;
        d.popoverTitle = popoverIntervalTitle(d,i);
        d.popoverContent = popoverIntervalContent(d,i);
      })
      .attr('x', function(d,i) { return scale(d.startPoint); })
      .attr('y', function(d,i) { return yScale(d.y) - 0.5 * margins.bar; })
      .attr('width', function(d,i) { return scale(d.endPoint) - scale(d.startPoint); })
      .attr('height', margins.bar);

    shapes
      .attr('x', function(d,i) { return scale(d.startPoint); })
      .attr('y', function(d,i) { return yScale(d.y) - 0.5 * margins.bar; })
      .attr('width', function(d,i) { return scale(d.endPoint) - scale(d.startPoint); })
      .attr('height', margins.bar)
      .style('fill', function(d,i) { return metadata[d.chromosome].color; })
      .style('stroke', function(d,i) { return d3.rgb(metadata[d.chromosome].color).darker(1); })
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlighted', false);
      })
      .on('mousemove', function(d,i) {
        var popover = d3.select('.popover');
        popover.select('.popover-title').html(d.popoverTitle);
        popover.select('.popover-content').html(d.popoverContent);
        popover.select('.popover-content span').style('color', d.color)
        popover
          .style('left', (d3.event.pageX - 0.91 *  popover.node().getBoundingClientRect().width / 2) + 'px')
          .style('top', (d3.event.pageY - popover.node().getBoundingClientRect().height - 3) + 'px')
          .classed('hidden', false)
          .style('display', 'block')
          .transition()
          .duration(5)
          .style('opacity', 1);
      });

    shapes.exit().remove();
  }

  function drawLocalConnections(container) {

    var connections = container.selectAll('path.connection').data(function(d,i) { return (localConnectionChromosomeBins[d.chromosome] || []); }, function(d,i) { return d.cid});

    connections.exit().remove();

    connections.attr('d', function(d,i) { return line(calculateConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); });

    connections
      .enter()
      .append('path')
      .attr('class', function(d,i) { return 'popovered connection local ' + d.type; })
      .style('clip-path','url(#clip)')
      .attr('d', function(d,i) { return line(calculateConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); })
      .each(function(d,i) {
        d.popoverTitle = popoverConnectionTitle(d,i);
        d.popoverContent = popoverConnectionContent(d,i);
      })
      .on('mouseover', function(d,i) {
        var connector = connectionBins[d.cid];
        d3.select(this).classed('highlighted', true);
        d3.selectAll('rect.shape').filter(function(e,j) { return [connector.source.iid,connector.sink.iid].includes(e.iid); }).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        var connector = connectionBins[d.cid];
        d3.select(this).classed('highlighted', false);
        d3.selectAll('rect.shape').filter(function(e,j) { return [connector.source.iid,connector.sink.iid].includes(e.iid); }).classed('highlighted', false);
      })
      .on('mousemove', function(d,i) {
        var popover = d3.select('.popover');
        popover.select('.popover-title').html(d.popoverTitle);
        popover.select('.popover-content').html(d.popoverContent);
        popover.select('.popover-content span').style('color', d.color)
        popover
          .style('left', (d3.event.pageX - 0.91 *  popover.node().getBoundingClientRect().width / 2) + 'px')
          .style('top', (d3.event.pageY - popover.node().getBoundingClientRect().height - 3) + 'px')
          .classed('hidden', false)
          .style('display', 'block')
          .transition()
          .duration(5)
          .style('opacity', 1);
      });
  }

  function drawInterChromosomeConnections(container) {

    var connections = container.selectAll('path.connection').data(function(d,i) { return (interChromosomeConnectionBins || []); }, function(d,i) { return d.cid});
 
    connections.exit().remove();

    connections.attr('d', function(d,i) { return line(calculateInterConnectorEndpoints(yScale, d, connectionBins[d.cid], panels)); });

    connections
      .enter()
      .append('path')
      .attr('class', function(d,i) { return 'popovered connection local ' + d.type; })
      .style('clip-path','url(#clipWidth)')
      .attr('d', function(d,i) { return line(calculateInterConnectorEndpoints(yScale, d, connectionBins[d.cid], panels)); })
      .each(function(d,i) {
        d.popoverTitle = popoverConnectionTitle(d,i);
        d.popoverContent = popoverConnectionContent(d,i);
      })
      .on('mouseover', function(d,i) {
        var connector = connectionBins[d.cid];
        d3.select(this).classed('highlighted', true);
        d3.selectAll('rect.shape').filter(function(e,j) { return [connector.source.iid,connector.sink.iid].includes(e.iid); }).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        var connector = connectionBins[d.cid];
        d3.select(this).classed('highlighted', false);
        d3.selectAll('rect.shape').filter(function(e,j) { return [connector.source.iid,connector.sink.iid].includes(e.iid); }).classed('highlighted', false);
      })
      .on('mousemove', function(d,i) {
        var popover = d3.select('.popover');
        popover.select('.popover-title').html(d.popoverTitle);
        popover.select('.popover-content').html(d.popoverContent);
        popover.select('.popover-content span').style('color', d.color)
        popover
          .style('left', (d3.event.pageX - 0.91 *  popover.node().getBoundingClientRect().width / 2) + 'px')
          .style('top', (d3.event.pageY - popover.node().getBoundingClientRect().height - 3) + 'px')
          .classed('hidden', false)
          .style('display', 'block')
          .transition()
          .duration(5)
          .style('opacity', 1);
      });
  }

  function drawLooseConnections(container) {

    var connections = container.selectAll('path.connection').data(function(d,i) { return (looseConnectionChromosomeBins[d.chromosome] || []); }, function(d,i) { return d.cid});

    connections.exit().remove();

    connections.attr('d', function(d,i) { return line(calculateLooseConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); });

    connections
      .enter()
      .append('path')
      .attr('class', function(d,i) { return 'popovered connection local ' + d.type; })
      .style('clip-path','url(#clip)')
      .attr('d', function(d,i) { return line(calculateLooseConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); })
      .each(function(d,i) {
        d.popoverTitle = popoverConnectionTitle(d,i);
        d.popoverContent = popoverConnectionContent(d,i);
      })
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlighted', false);
      })
      .on('mousemove', function(d,i) {
        var popover = d3.select('.popover');
        popover.select('.popover-title').html(d.popoverTitle);
        popover.select('.popover-content').html(d.popoverContent);
        popover.select('.popover-content span').style('color', d.color)
        popover
          .style('left', (d3.event.pageX - 0.91 *  popover.node().getBoundingClientRect().width / 2) + 'px')
          .style('top', (d3.event.pageY - popover.node().getBoundingClientRect().height - 3) + 'px')
          .classed('hidden', false)
          .style('display', 'block')
          .transition()
          .duration(5)
          .style('opacity', 1);
      });
  }

  function drawLocalInterConnections(container) {

    var connections = container.selectAll('path.connection').data(function(d,i) { return (localInterChromosomeConnectionBins[d.chromosome] || []); }, function(d,i) { return d.cid});

    connections.exit().remove();

    connections.attr('d', function(d,i) { return line(calculateLocalInterConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); });

    connections
      .enter()
      .append('path')
      .attr('class', function(d,i) { return 'popovered connection local ' + d.type; })
      .style('clip-path','url(#clip)')
      .attr('d', function(d,i) { return line(calculateLocalInterConnectorEndpoints(yScale, d, connectionBins[d.cid], d3.select(this.parentNode).datum())); })
      .each(function(d,i) {
        d.popoverTitle = popoverConnectionTitle(d,i);
        d.popoverContent = popoverConnectionContent(d,i);
        var chromosomeObject = d3.select(this.parentNode).datum();
        d.outerChromosome = (connectionBins[d.cid].source.chromosome === chromosomeObject.chromosome) ? connectionBins[d.cid].sink.chromosome : connectionBins[d.cid].source.chromosome;
      })
      .style('stroke', function(d,i) { return d3.rgb(metadata[d.outerChromosome].color).darker(1); })
      .on('click', function(d,i) {
        /*
        refreshPanel(0, connectionBins[d.cid].source.chromosome);
        refreshPanel(1, d.outerChromosome);
        refreshPanel(2, "17");
        */
        var column = 1;
        d3.select('.control-container.column-' + column).selectAll('circle').classed('selected', false);
        d3.select('.control-container.column-' + column).selectAll('circle').filter(function(e,j) { return e.chromoObject.chromosome === d.outerChromosome }).classed('selected', true);
        panels[column] = Object.assign({}, metadata[d.outerChromosome], {column: column});
        interChromosomeConnectionBins = getInterChromosomeConnectionBins(data, panels, connectionBins);
        localInterChromosomeConnectionBins = getLocalInterChromosomeConnectionBins(data, panels, connectionBins);
        updatePanels(panels);
        updateLegend(panels);
      })
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlighted', false);
      })
      .on('mousemove', function(d,i) {
        var popover = d3.select('.popover');
        popover.select('.popover-title').html(d.popoverTitle);
        popover.select('.popover-content').html(d.popoverContent);
        popover.select('.popover-content span').style('color', d.color)
        popover
          .style('left', (d3.event.pageX - 0.91 *  popover.node().getBoundingClientRect().width / 2) + 'px')
          .style('top', (d3.event.pageY - popover.node().getBoundingClientRect().height - 3) + 'px')
          .classed('hidden', false)
          .style('display', 'block')
          .transition()
          .duration(5)
          .style('opacity', 1);
      });
  }

  // Callback when brushing is finished
  function brushed() {
    if (d3.event.sourceEvent && d3.event.sourceEvent.type === 'zoom') return; // ignore brush-by-zoom
    var s = d3.event.selection || [0, panelContainerWidth];
    var brushData = d3.select(this).datum();
    var domain = s.map(brushData.scale2.invert, brushData.scale2);
    brushData.scale.domain(domain);
    brushData.panelScale.domain(domain);
    var panel = d3.selectAll('.panel-container').filter(function(d,i) { return d.column === brushData.column && d.chromosome === brushData.chromosome });
    panel.select('.axis--x').call(brushData.axis).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start');
    panel.select('.zoom').call(brushData.zoom.transform, d3.zoomIdentity.scale(panelContainerWidth / (s[1] - s[0])).translate(-s[0], 0));
    var intervals = data.intervals.filter(function(d,i) { return (d.chromosome === brushData.chromosome)});
    drawIntervals(panel.select('g.shapes-container'), brushData.scale, intervals);
    drawLocalConnections(panel.select('g.local-connections-container'));
    drawLooseConnections(panel.select('g.loose-connections-container'));
    drawLocalInterConnections(panel.select('g.local-inter-connections-container'));
    drawInterChromosomeConnections(svg.select('g.inter-chromosome-connections-container'));
  }

  // Callback when the panel is zoomed
  function zoomed(panelData) {
    if (d3.event.sourceEvent && d3.event.sourceEvent.type === 'brush') return; // ignore zoom-by-brush
    var t = d3.event.transform;
    var panel = d3.selectAll('.panel-container').filter(function(d,i) { return d.column === panelData.column && d.chromosome === panelData.chromosome });
    var chromo = d3.selectAll('.chromosome-container').filter(function(d,i) { return d.column === panelData.column && d.chromosome === panelData.chromosome });
    var domain = t.rescaleX(panelData.scale2).domain();
    panelData.scale.domain(domain);
    panelData.panelScale.domain(domain);
    panel.select('.axis--x').call(panelData.axis).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start');
    chromo.select('.brush').call(panelData.brush.move, panelData.scale.range().map(t.invertX, t));
    var intervals = data.intervals.filter(function(d,i) { return (d.chromosome === panelData.chromosome)});
    drawIntervals(panel.select('g.shapes-container'), panelData.scale, intervals);
    drawLocalConnections(panel.select('g.local-connections-container'));
    drawLooseConnections(panel.select('g.loose-connections-container'));
    drawLocalInterConnections(panel.select('g.local-inter-connections-container'));
    drawInterChromosomeConnections(svg.select('g.inter-chromosome-connections-container'));
  }

  function refreshPanel(theColumn, theChromosome) {
    d3.select('.control-container.column-' + theColumn).selectAll('circle').classed('selected', function(d,i) { return d.chromoObject.chromosome === theChromosome });
    panels[theColumn] = Object.assign({}, metadata[theChromosome], {column: theColumn});
    interChromosomeConnectionBins = getInterChromosomeConnectionBins(data, panels, connectionBins);
    localInterChromosomeConnectionBins = getLocalInterChromosomeConnectionBins(data, panels, connectionBins);
    updatePanels(panels);
    updateLegend(panels);
  }
}

// Act upon window resize
function throttle() {
  window.clearTimeout(throttleTimer);
  throttleTimer = window.setTimeout(function() {
    draw();
  }, 200);
}

// Remove any other open popovers
$(document).on('mousemove', function(event) {
  if (!$(event.target).is('.popovered')) {
    d3.select('.popover').transition().duration(5)
      .style('opacity', 0);
  }
});

