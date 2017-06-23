class BrushContainer {

  constructor(frame) {
    this.frame = frame;
    this.reset();
  }

  reset() {
    this.activeId = null;
    this.originalSelection = null;
    this.currentSelection = null;
    this.totalBrushWidth = this.frame.genomeScale.range()[1] - this.frame.genomeScale.range()[0];
    this.otherSelections = [];
    this.fragments = [];
    this.visibleFragments = [];
    this.visibleIntervals = [];
    this.panelWidth = 0;
    this.panelHeight = 0;
  }

  render() {
    this.reset();
    this.createBrush();
    this.update();
  }

  deleteBrush() {
    this.fragments = this.fragments.filter(fragment => fragment.id !== this.activeId);
    this.update();
  }

  createBrush() {
    var self = this;
    var brush = d3.brushX()
      .extent([[0, 0], [this.totalBrushWidth, this.frame.margins.brushes.height]])
      .on('start', function() {
        // brush starts here
        self.originalSelection = d3.event.selection;
      })
      .on('brush', function() {
        // brushing happens here

        // ignore brush-by-zoom
        if (d3.event.sourceEvent && d3.event.sourceEvent.type === 'zoom') return;

        // Only transition after input.
        if (!d3.event || !d3.event.sourceEvent || (d3.event.sourceEvent.type === 'brush')) return;

        let fragment = d3.select(this).datum();
        self.activeId = d3.select(this).datum().id;
        let originalSelection = fragment.selection;
        let currentSelection = d3.event.selection;
        let selection = Object.assign([], currentSelection);
        let node;

        // read the current state of all the self.fragments before you start checking on collisions
        self.otherSelections = self.fragments.filter((d, i) => (d.selection !== null) && (d.id !== self.activeId)).map((d, i) => {
          node = d3.select('#brush-' + d.id).node();
          return node && d3.brushSelection(node); 
        });

        // calculate the lower allowed selection edge this brush can move
        let lowerEdge = d3.max(self.otherSelections.filter((d, i) => (d.selection !== null))
          .filter((d, i) => originalSelection && (d[0] <= originalSelection[0]) && (originalSelection[0] <= d[1]))
          .map((d, i) => d[1]));

        // calculate the upper allowed selection edge this brush can move
        let upperEdge = d3.min(self.otherSelections.filter((d, i) => (d.selection !== null))
          .filter((d, i) => originalSelection && (d[1] >= originalSelection[0]) && (originalSelection[1] <= d[1]))
          .map((d, i) => d[0]));

        // if there is an upper edge, then set this to be the upper bound of the current selection
        if ((upperEdge !== undefined) && (selection[1] >= upperEdge)) {
          selection[1] = upperEdge;
          selection[0] = d3.min([selection[0], upperEdge - 1]);
        } 

        // if there is a lower edge, then set this to the be the lower bound of the current selection
        if ((lowerEdge !== undefined) && (selection[0] <= lowerEdge)) {
          selection[0] = lowerEdge;
          selection[1] = d3.max([selection[1], lowerEdge + 1]);
        }

        // move the brush to stay within the allowed bounded selection zone
        if ((selection !== undefined) && (selection !== null) && (selection[1] !== selection[0])) {
          d3.select(this).call(fragment.brush.move, selection);
        }

        // finally, update the chart with the selection in question
        self.update();
      })
      .on('end', function() {
        // ignore brush-by-zoom
        if (d3.event.sourceEvent && d3.event.sourceEvent.type === 'zoom') return;
        
        // Only transition after input.
        if (!d3.event.sourceEvent) return;

        // Ignore empty selections.
        if (!d3.event.selection) return;

        // Figure out if our latest brush has a selection
        let lastBrushID = self.fragments[self.fragments.length - 1].id;
        let lastBrush = d3.select('#brush-' + lastBrushID).node();
        let selection = d3.brushSelection(lastBrush);

        // If it does, that means we need another one
        if (selection && selection[0] !== selection[1]) {
          self.createBrush();
        }

        // finally, update the chart with the selection in question
        self.update();
    });

    this.fragments.push(new Fragment(brush));
  }

  update() {

    // first recalculate the current selections
    this.updateFragments();

    // Draw the notes of the fragments
    this.renderFragmentsNote();

    // update clipPath
    this.renderClipPath();

    // draw the brushes
    this.renderBrushes();

    // Draw the panel rectangles
    this.renderPanels();

    // Draw the chromosome axis
    this.renderChromoAxis();

    // Draw the intervals
    this.renderIntervals();

    // Draw the interconnections
    this.renderInterconnections();
  }

  updateFragments() {
    let node;
    this.visibleFragments = [];
    this.visibleIntervals = [];
    this.connections = [];

    let frameConnections = Object.assign([], this.frame.connections);

  // delete any brushes that have a zero selection size
    this.fragments = this.fragments.filter((d, i) => (d.selection === null) || (d.selection[0] !== d.selection[1]));

  // filter the brushes that are visible on the screen
    this.fragments.forEach((fragment, i) => {
      node = d3.select('#brush-' + fragment.id).node();
      fragment.selection = node && d3.brushSelection(node);
      fragment.domain = fragment.selection && fragment.selection.map(this.frame.genomeScale.invert,this.frame.genomeScale);
      if (fragment.selection) {
        this.visibleFragments.push(Object.assign({}, fragment));
      }
    });

    // determine the new Panel Width
    this.panelWidth = (this.frame.width - (this.visibleFragments.length - 1) * this.frame.margins.panels.gap) / this.visibleFragments.length;
    this.panelHeight = this.frame.height - this.frame.margins.panels.upperGap + this.frame.margins.top;

    // now sort the visible self.fragments from smallest to highest
    this.visibleFragments = Object.assign([], this.visibleFragments.sort((x, y) => d3.ascending(x.selection[0], y.selection[0])));

    // Determine the panel parameters for rendering
    this.visibleFragments.forEach((d, i) => {
      d.panelWidth = this.panelWidth;
      d.panelHeight = this.panelHeight;
      d.range = [i * (d.panelWidth + this.frame.margins.panels.gap), (i + 1) * d.panelWidth + i * this.frame.margins.panels.gap];
      d.scale = d3.scaleLinear().domain(d.domain).range(d.range);
      d.innerScale = d3.scaleLinear().domain(d.domain).range([0, d.panelWidth]);
      d.zoom = d3.zoom().scaleExtent([1, Infinity]).translateExtent([[0, 0], [this.frame.width, d.panelHeight]]).extent([[0, 0], [this.frame.width, d.panelHeight]]).on('zoom', () => this.zoomed(d));
      d.chromoAxis = Object.keys(this.frame.chromoBins)
        .map(x => this.frame.chromoBins[x])
        .filter(chromo => chromo.contains(d.domain))
        .map((chromo, j) => {
          let domainStart = ((d.domain[0] < chromo.scaleToGenome.range()[1]) && (d.domain[0] >= chromo.scaleToGenome.range()[0])) ? chromo.scaleToGenome.invert(d.domain[0]) : chromo.scaleToGenome.domain()[0];
          let domainEnd   = ((d.domain[1] < chromo.scaleToGenome.range()[1]) && (d.domain[1] >= chromo.scaleToGenome.range()[0])) ? chromo.scaleToGenome.invert(d.domain[1]) : chromo.scaleToGenome.domain()[1];
          let rangeWidth = d.innerScale(chromo.scaleToGenome(domainEnd)) - d.innerScale(chromo.scaleToGenome(domainStart));
          let scale = d3.scaleLinear().domain([domainStart, domainEnd]).range([0, rangeWidth]);
          let axisBottom = d3.axisBottom(scale).ticks(d3.max([d3.min([Math.round(rangeWidth / 25), 10]),1]), 's');
          return {identifier: Misc.guid, transform: 'translate(' + [d.innerScale(chromo.scaleToGenome(domainStart)), 0] + ')',
            labelTopTranslate: 'translate(' + [0.5 * (d.innerScale(chromo.scaleToGenome(domainEnd)) - d.innerScale(chromo.scaleToGenome(domainStart))), - this.frame.margins.panels.label] + ')',
            chromo: chromo, scale: scale, rangeWidth: rangeWidth, axisBottom: axisBottom};
      });
      // filter the intervals
      d.visibleIntervals = [];
      this.frame.intervals
      .filter((e, j) => ((e.startPlace <= d.domain[1]) && (e.startPlace >= d.domain[0])) || ((e.endPlace <= d.domain[1]) && (e.endPlace >= d.domain[0])) 
        || (((d.domain[1] <= e.endPlace) && (d.domain[1] >= e.startPlace)) || ((d.domain[0] <= e.endPlace) && (d.domain[0] >= e.startPlace))))
      .forEach((interval, j) => {
        interval.identifier = Misc.guid;
        interval.range = [d3.max([0, d.innerScale(interval.startPlace)]), d.innerScale(interval.endPlace)];
        interval.shapeWidth = interval.range[1] - interval.range[0];
        interval.fragment = d;
        d.visibleIntervals.push(interval);
      });
      // filter the connections on same fragment
      frameConnections
        .filter((e, j) => (!e.source || ((e.source.place <= d.domain[1]) && (e.source.place >= d.domain[0]))) && (!e.sink || ((e.sink.place <= d.domain[1]) && (e.sink.place >= d.domain[0]))))
        .forEach((connection, j) => {
          if (connection.source) {
            connection.source.scale = d.scale;
            connection.source.fragment = d;
          }
          if (connection.sink) {
            connection.sink.scale = d.scale;
            connection.sink.fragment = d;
          }
          connection.touchScale = d.scale;
          connection.identifier = Misc.guid;
          this.connections.push(connection);
        });
    });
    // filter the connections between the visible fragments
    k_combinations(this.visibleFragments, 2).forEach((pair, i) => {
      frameConnections
        .filter((e, j) => (e.type !== 'LOOSE') 
          && (((e.source.place <= pair[0].domain[1]) && (e.source.place >= pair[0].domain[0]) && (e.sink.place <= pair[1].domain[1]) && (e.sink.place >= pair[1].domain[0]))
          ||((e.source.place <= pair[1].domain[1]) && (e.source.place >= pair[1].domain[0]) && (e.sink.place <= pair[0].domain[1]) && (e.sink.place >= pair[0].domain[0]))))
        .forEach((connection, j) => {
          if ((connection.source.place <= pair[0].domain[1]) && (connection.source.place >= pair[0].domain[0])) {
            connection.source.scale = pair[0].scale;
            connection.source.fragment = pair[0];
          } else {
            connection.source.scale = pair[1].scale;
            connection.source.fragment = pair[1];
          }
          if ((connection.sink.place <= pair[0].domain[1]) && (connection.sink.place >= pair[0].domain[0])) {
            connection.sink.scale = pair[0].scale;
            connection.sink.fragment = pair[0];
          } else {
            connection.sink.scale = pair[1].scale;
            connection.sink.fragment = pair[1];
          }
          connection.identifier = Misc.guid;
          this.connections.push(connection);
        });
    });
    let visibleConnections = Object.assign([], this.connections).map((d, i) => d.cid);
    this.visibleFragments.forEach((fragment, i) => {
      frameConnections
        .filter((e, j) => { return (e.type !== 'LOOSE') && (!visibleConnections.includes(e.cid))
          && (((e.source.place <= fragment.domain[1]) && (e.source.place >= fragment.domain[0]))
          ||((e.sink.place <= fragment.domain[1]) && (e.sink.place >= fragment.domain[0])))})
        .forEach((con, j) => {
          let connection = Object.assign(new Connection(con), con);
          connection.locateAnchor(fragment);
          this.connections.push(connection);
        });
    });
  }

  zoomed(fragment) {
    var self = this;
    if (d3.event.sourceEvent && d3.event.sourceEvent.type === 'brush') return; // ignore zoom-by-brush
    // set this brush as active

    // Get the generated domain upon zoom 
    let t = d3.event.transform;
    let zoomedDomain = t.rescaleX(this.frame.genomeScale).domain();
    let domain = Object.assign([], zoomedDomain);
    //console.log(t, t.rescaleX(this.frame.genomeScale).domain(), t.rescaleX(fragment.innerScale).domain())

    // Calculate the other domains and the domain bounds for the current brush
    let otherDomains = this.fragments.filter((d, i) => (d.selection !== null) && (d.id !== fragment.id)).map((d, i) => d.domain);
    let lowerBound = d3.max(otherDomains.filter((d, i) => fragment.domain && (d[1] <= fragment.domain[0])).map((d, i) => d[1]));
    let upperBound = d3.min(otherDomains.filter((d, i) => fragment.domain && (d[0] >= fragment.domain[1])).map((d, i) => d[0]));

    // if there is an upper bound, set this to the maximum allowed limit
    if ((upperBound !== undefined) && (domain[1] >= upperBound)) {
      domain[1] = upperBound;
      domain[0] = d3.min([domain[0], upperBound - 1]);
    } 
    // if there is a lower bound, set this to the lowest allowed limit
    if ((lowerBound !== undefined) && (domain[0] <= lowerBound)) {
      domain[0] = lowerBound;
      domain[1] = d3.max([domain[1], lowerBound + 1]);
    }

    // update the current brush
    fragment.scale.domain(domain);
    let selection = [this.frame.genomeScale(domain[0]), this.frame.genomeScale(domain[1])];
    d3.select('#brush-' + fragment.id).call(fragment.brush.move, selection);

    // update the data
    this.updateFragments();

    // update the interconnections
    this.renderInterconnections();
    
    //update the panel axis Top
    this.frame.panelsChromoAxisContainerTop.selectAll('g.axis')
      .data(this.visibleFragments,  (d, i) => d.id)
      .each(function(d,i) {
        d3.select(this).select('rect').attr('width', (e, j) => d.rangeWidth);
        d3.select(this).select('text.label-legend').attr('transform', (e, j) => d.labelTopTranslate); 
      });

    //update the panel axis Top
    this.frame.panelsChromoAxisContainerBottom.selectAll('g.axis')
      .data(this.visibleFragments,  (d, i) => d.id)
      .each(function(d,i) { 
        d3.select(this).call(d.axisBottom).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start'); 
      });

    // update the chromosome axis
    this.renderChromoAxis();

    // update the intervals
    this.renderIntervals();
  }

  renderClipPath() {
    if (this.visibleFragments.length > 0) {
      this.frame.svgFilter.renderClipPath(this.panelWidth + this.frame.margins.panels.widthOffset, this.panelHeight);
    }
  }

  renderBrushes() {
    var self = this;

    let brushSelection = this.frame.brushesContainer.selectAll('.brush')
      .data(this.fragments,  (d, i) => d.id);

    // Set up new brushes
    brushSelection
      .enter()
      .insert('g', '.brush')
      .attr('class', 'brush')
      .attr('id', (d, i) => 'brush-' + d.id)
      .each(function(fragment) {
        //call the brush
        d3.select(this).call(fragment.brush);
      });

    // update the brushes
    brushSelection
      .each(function (fragment){
        d3.select(this)
          .attr('class', 'brush')
          .classed('highlighted', (d, i) => d.id === self.activeId)
          .selectAll('.overlay')
          .style('pointer-events',(d, i) => {
            let brush = fragment.brush;
            if (fragment.id === self.fragments[self.fragments.length - 1].id && brush !== undefined) {
              return 'all';
            } else {
              return 'none';
            }
          });
      })

    // exit the brushes
    brushSelection
      .exit()
      .remove();
  }

  renderPanels() {
    let self = this;
    let correctionOffset = 1; // used for aligning the rectenges on the y Axis lines

    // Draw the panel rectangles
    let panelRectangles = this.frame.panelsContainer.selectAll('rect.panel')
      .data(this.visibleFragments,  (d, i) => d.id);

    panelRectangles
      .enter()
      .append('rect')
      .attr('class', 'panel')
      .style('clip-path','url(#clip)')
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')')
      .attr('width', (d, i) => d.panelWidth + this.frame.margins.panels.widthOffset)
      .attr('height', (d, i) => d.panelHeight + correctionOffset)
      .each(function(d,i) {
        d3.select(this)
          .call(d.zoom.transform, d3.zoomIdentity
          .scale(self.frame.width / (d.selection[1] - d.selection[0]))
          .translate(-d.selection[0], 0));
      })
      .on('click', (d,i) => {
        this.activeId = d.id;
        this.frame.brushesContainer.selectAll('.brush').classed('highlighted', false);
        d3.select('#brush-' + d.id).classed('highlighted', true);
      })

    panelRectangles
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')')
      .attr('width', (d, i) => d.panelWidth + this.frame.margins.panels.widthOffset)
      .each(function(d,i) {
        d3.select(this).call(d.zoom)
         .call(d.zoom.transform, d3.zoomIdentity
         .scale(self.frame.width / (d.selection[1] - d.selection[0]))
         .translate(-d.selection[0], 0));
      });

    panelRectangles
      .exit()
      .remove();
  }

  renderChromoAxis() {
    let self = this;

    //Chromo Axis Top
    let chromoAxisContainer = this.frame.panelsChromoAxisContainerTop.selectAll('g.chromo-axis-container')
      .data(this.visibleFragments,  (d, i) => d.id);

    chromoAxisContainer
      .enter()
      .append('g')
      .attr('class', 'chromo-axis-container')
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    chromoAxisContainer
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    chromoAxisContainer
      .exit()
      .remove();

    let chromoAxis = chromoAxisContainer.selectAll('g.chromo-legend')
      .data((d, i) => d.chromoAxis, (d, i) => d.identifier);

    chromoAxis
      .enter()
      .append('g')
      .attr('class', 'chromo-legend')
      .attr('transform', (d, i) => d.transform)
      .each(function(d,i) {
        d3.select(this).append('rect').attr('width', (e, j) => d.rangeWidth).attr('y', -self.frame.margins.panels.legend).attr('height', self.frame.margins.panels.legend).style('fill', (e, j) => "url('#gradient" + d.chromo.chromosome +"')");
        d3.select(this).append('text').attr('class', 'label-chromosome').attr('transform', (e, j) => d.labelTopTranslate).text((e, j) => d.chromo.chromosome);
      })

    chromoAxis
      .attr('transform', (d, i) => d.transform)
      .each(function(d,i) { 
        d3.select(this).select('rect').attr('width', (e, j) => d.rangeWidth);
        d3.select(this).select('text.label-chromosome').attr('transform', (e, j) => d.labelTopTranslate);
      });

    chromoAxis
      .exit()
      .remove();

    // Chromo Axis Bottom
    let chromoAxisContainerBottom = this.frame.panelsChromoAxisContainerBottom.selectAll('g.chromo-axis-container-bottom')
      .data(this.visibleFragments,  (d, i) => d.id);

    chromoAxisContainerBottom
      .enter()
      .append('g')
      .attr('class', 'chromo-axis-container-bottom')
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    chromoAxisContainerBottom
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    chromoAxisContainerBottom
      .exit()
      .remove();

    let chromoAxisBottom = chromoAxisContainerBottom.selectAll('g.chromo-axis-bottom')
      .data((d, i) => d.chromoAxis, (d, i) => d.identifier);

    chromoAxisBottom
      .enter()
      .append('g')
      .attr('class', 'chromo-axis-bottom axis axis--x')
      .attr('transform', (d, i) => d.transform)
      .each(function(d,i) {
        d3.select(this).call(d.axisBottom).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start').style('fill', (e, j) => d.chromo.color);
      });

    chromoAxisBottom
      .attr('transform', (d, i) => d.transform)
      .each(function(d,i) { 
        d3.select(this).call(d.axisBottom).selectAll('text').attr('transform', 'rotate(45)').style('text-anchor', 'start').style('fill', (e, j) => d.chromo.color);
      });

    chromoAxisBottom
      .exit()
      .remove();
  }

  renderIntervals() {
    // create the g elements containing the intervals
    let shapesPanels = this.frame.shapesContainer.selectAll('g.shapes-panel')
      .data(this.visibleFragments,  (d, i) => d.id);

    shapesPanels
      .enter()
      .append('g')
      .attr('class', 'shapes-panel')
      .style('clip-path','url(#clip)')
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    shapesPanels
      .attr('transform', (d, i) => 'translate(' + [d.range[0], 0] + ')');

    shapesPanels
      .exit()
      .remove();

    // add the actual intervals as rectangles
    let shapes = shapesPanels.selectAll('rect.shape')
      .data((d, i) => d.visibleIntervals, (d, i) =>  d.identifier);

    shapes
      .enter()
      .append('rect')
      .attr('id', (d, i) => d.identifier)
      .attr('class', 'popovered shape')
      .attr('transform', (d, i) => 'translate(' + [d.range[0], this.frame.yScale(d.y) - 0.5 * this.frame.margins.intervals.bar] + ')')
      .attr('width', (d, i) => d.shapeWidth)
      .attr('height', this.frame.margins.intervals.bar)
      .style('fill', (d, i) => d.color)
      .style('stroke', (d, i) => d3.rgb(d.color).darker(1))
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlighted', false);
      })
      .on('mousemove', (d,i) => this.loadPopover(d))
      .on('dblclick', (d,i) => {
        let fragment = d.fragment;
        let lambda = (fragment.panelWidth - 2 * this.frame.margins.intervals.gap) / (d.endPlace - d.startPlace);
        let domainOffset = this.frame.margins.intervals.gap / lambda;
        fragment.domain = [d.startPlace - domainOffset, d.endPlace + domainOffset];
        fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
        d3.select('#brush-' + fragment.id).call(fragment.brush.move, fragment.selection);
        this.update();
      });

    shapes
      .attr('id', (d, i) => d.identifier)
      .attr('transform', (d, i) => 'translate(' + [d.range[0], this.frame.yScale(d.y) - 0.5 * this.frame.margins.intervals.bar] + ')')
      .attr('width', (d, i) => d.shapeWidth)
      .style('fill', (d, i) => d.color)
      .style('stroke', (d, i) => d3.rgb(d.color).darker(1));

    shapes
      .exit()
      .remove();
  }

  renderInterconnections() {

    let connections = this.frame.connectionsContainer.selectAll('path.connection')
      .data(this.connections, (d,i) => d.identifier);
 
    connections.exit().remove();

    connections
      .attr('class', (d,i) => d.styleClass)
      .style('fill', (d, i) => d.fill)
      .style('stroke', (d, i) => d.stroke)
      .attr('transform', (d,i) => d.transform)
      .attr('d', (d,i) => d.render);

    connections
      .enter()
      .append('path')
      .attr('id', (d,i) => d.identifier)
      .attr('class', (d,i) => d.styleClass)
      .attr('transform', (d,i) => d.transform)
      .style('fill', (d, i) => d.fill)
      .style('stroke', (d, i) => d.stroke)
      .attr('d', (d,i) =>  d.render)
      .on('mouseover', function(d,i) {
        d3.select(this).classed('highlighted', true);
      })
      .on('mouseout', function(d,i) {
        d3.select(this).classed('highlighted', false);
      })
      .on('mousemove', (d,i) => this.loadPopover(d))
      .on('dblclick', (d,i) => {
        if (d.kind === 'ANCHOR') {
          this.createBrush();
          let fragment = this.fragments[this.fragments.length - 1];
          fragment.domain = [0.99 * d.otherEnd.interval.startPlace, 1.01 * d.otherEnd.interval.endPlace];
          fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
          this.update();
          fragment = d3.select('#brush-' + fragment.id).datum();
          let lambda = (this.panelWidth - 2 * this.frame.margins.intervals.gap) / (d.otherEnd.interval.endPlace - d.otherEnd.interval.startPlace);
          let domainOffset = this.frame.margins.intervals.gap / lambda;
          fragment.domain = [d.otherEnd.interval.startPlace - domainOffset, d.otherEnd.interval.endPlace + domainOffset];
          fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
          d3.select('#brush-' + fragment.id).call(fragment.brush.move, fragment.selection);
          this.update();
        } else {
          if (d.source.fragment.id === d.sink.fragment.id) {
            let fragment = d.source.fragment;
            let lambda = (fragment.panelWidth - 2 * this.frame.margins.intervals.gap) / Math.abs(d.source.place - d.sink.place);
            let domainOffset = this.frame.margins.intervals.gap / lambda;
            fragment.domain = [d3.min([d.source.place, d.sink.place]) - domainOffset, d3.max([d.source.place, d.sink.place]) + domainOffset];
            fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
            d3.select('#brush-' + fragment.id).call(fragment.brush.move, fragment.selection);
            this.update();
          } else {
            // first align the source interval
            let fragment = d.source.fragment;
            let lambda = (fragment.panelWidth - 2 * this.frame.margins.intervals.gap) / (d.source.interval.endPlace - d.source.interval.startPlace);
            let domainOffset = this.frame.margins.intervals.gap / lambda;
            fragment.domain = [d.source.interval.startPlace - domainOffset, d.source.interval.endPlace + domainOffset];
            fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
            d3.select('#brush-' + fragment.id).call(fragment.brush.move, fragment.selection);
            this.update();
            // second align the sink interval
            fragment = d.sink.fragment;
            lambda = (fragment.panelWidth - 2 * this.frame.margins.intervals.gap) / (d.sink.interval.endPlace - d.sink.interval.startPlace);
            domainOffset = this.frame.margins.intervals.gap / lambda;
            fragment.domain = [d.sink.interval.startPlace - domainOffset, d.sink.interval.endPlace + domainOffset];
            fragment.selection = [this.frame.genomeScale(fragment.domain[0]), this.frame.genomeScale(fragment.domain[1])];
            d3.select('#brush-' + fragment.id).call(fragment.brush.move, fragment.selection);
            this.update();
          }
        }
      });
  }

  renderFragmentsNote() {
    let note = this.visibleFragments.map((d, i) => d.chromoAxis.map((e, j) => {
      return (e.chromo.chromosome + ':' + Math.floor(e.scale.domain()[0]) + '-' + Math.floor(e.scale.domain()[1]));
    }).join(' ')).join(' | ');
    d3.select('#fragmentsNote').text(note);
  }

  loadPopover(d) {
    var popover = d3.select('.popover');
    popover.select('.popover-title').html(d.popoverTitle);
    popover.select('.popover-content').html(d.popoverContent);
    popover.select('.popover-content span').style('color', d.color)
    popover
      .style('left', (d3.event.pageX - 1.0 *  popover.node().getBoundingClientRect().width / 2) + 'px')
      .style('top', (d3.event.pageY - 1.31 * popover.node().getBoundingClientRect().height - 3) + 'px')
      .classed('hidden', false)
      .style('display', 'block')
      .transition()
      .duration(5)
      .style('opacity', 1);
  }
  
}