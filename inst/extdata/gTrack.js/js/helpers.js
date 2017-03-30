// The hashmap mapping the chromosome id to the chromosome metadata
function getMetadata(dataArray) {
  return dataArray.metadata.reduce(function(hash, elem) { hash[elem.chromosome] = elem; return hash }, {});
}

// The hashmap mapping each interval id to the interval itself
function getIntervalBins(dataArray) {
  return dataArray.intervals.reduce(function(hash, elem) { hash[elem.iid] = elem; return hash }, {});
}
// merge and deduplicate array in the form [ [1, 2, 3], [101, 2, 1, 10], [2, 1] ]
function getArraysIntersection(array1, array2) {
  array1.filter(function(n) {
    return array2.indexOf(n) !== -1;
  })
}
// The hashmap mapping the connection id to the connection and its intervals
function getConnectionBins(dataArray, intervalBins) {
  return dataArray.connections.reduce(function(hash, elem) { 
    hash[elem.cid] = {
      connection: elem,
      source: intervalBins[elem.source && Math.abs(elem.source)],
      sink: intervalBins[elem.sink && Math.abs(elem.sink)]
    }; 
    return hash; }, {});
}
// The hashmap mapping the connection id to the connections within the same chromosome
function getLocalConnectionBins(dataArray, connectionBins) {
  return dataArray.connections.reduce(function(hash, elem) {
  var source = connectionBins[elem.cid].source;
  var sink = connectionBins[elem.cid].sink;
  if ((elem.type !== 'LOOSE') && (source.chromosome === sink.chromosome)) {
    if (!hash[source.chromosome]) {
      hash[source.chromosome] = [];
    }
    hash[source.chromosome].push(elem);
  }
  return hash; }, {});
}
// The hashmap mapping the connection id to the connections between intervals in different chromosomes selected on the panels
function getInterChromosomeConnectionBins(dataArray, panels, connectionBins) {
  var panelChromosomes = panels.map(function(d,i) { return d.chromosome; });
  var source, sink, verdict;
  return dataArray.connections.filter(function(elem, index) {
    source = connectionBins[elem.cid].source;
    sink = connectionBins[elem.cid].sink;
    verdict = ((elem.type !== 'LOOSE') && (source.chromosome !== sink.chromosome) && (panelChromosomes.includes(source.chromosome)) && (panelChromosomes.includes(sink.chromosome)));
    //verdict = verdict && ([source.chromosome, sink.chromosome].filter(function(n) {return [panelChromosomes[0], panelChromosomes[panelChromosomes.length - 1]].indexOf(n) !== -1;}).length < 2);
    return verdict;
  });
}
// The hashmap mapping the connection id to the connections between intervals in different chromosomes, not selected on the panels
function getLocalInterChromosomeConnectionBins(dataArray, panels, connectionBins) {
  var panelChromosomes = panels.map(function(d,i) { return d.chromosome; });
  return dataArray.connections.reduce(function(hash, elem) {
  var source = connectionBins[elem.cid].source;
  var sink = connectionBins[elem.cid].sink;
  if ((elem.type !== 'LOOSE') && (source.chromosome !== sink.chromosome) && !((panelChromosomes.includes(source.chromosome)) && (panelChromosomes.includes(sink.chromosome)))) {
    if (!hash[source.chromosome]) {
      hash[source.chromosome] = [];
    }
    if (!hash[sink.chromosome]) {
      hash[sink.chromosome] = [];
    }
    hash[source.chromosome].push(elem);
    hash[sink.chromosome].push(elem);
  }
  return hash; }, {});
}
// The hashmap mapping the connection id to the loose connections within the same chromosome
function getLooseConnectionBins(dataArray, connectionBins) {
  return dataArray.connections.reduce(function(hash, elem) {
  var source = connectionBins[elem.cid].source;
  var sink = connectionBins[elem.cid].sink;
  if (!source && sink) {
    if (!hash[sink.chromosome]) {
      hash[sink.chromosome] = [];
    }
    hash[sink.chromosome].push(elem);
  }
  if (source && !sink) {
    if (!hash[source.chromosome]) {
      hash[source.chromosome] = [];
    }
    hash[source.chromosome].push(elem);
  }
  return hash; }, {});
}
// The title for the popover on the intervals
function popoverIntervalTitle(d,i) {
  return 'Interval #' + d.title;
}

// The content for the popover of the intervals
function popoverIntervalContent(d,i) {
  var content = '';
  [{label: 'Chromosome', value: d.chromosome}, {label: 'Jabba', value: d.y}, {label: 'Start Point', value: d3.format(',')(d.startPoint)},
   {label: 'End Point', value: d3.format(',')(d.endPoint)}, {label: 'Interval Length', value: d3.format(',')(d.intervalLength)}, {label: 'Strand', value: d.strand}].forEach(function(e,j) {
     content += '<tr><td class="table-label" align="left" width="200" valign="top"><strong>' + e.label + ':</strong></td><td class="table-value" width="100" align="right" valign="top">' + e.value + '</td></tr>';
   });
  return '<div class="row"><div class="col-lg-12"><table width="0" border="0" align="left" cellpadding="0" cellspacing="0"><tbody>' + content + '</tbody></table></div></div>';
}

// The title for the popover on the connections
function popoverConnectionTitle(d,i) {
  return 'Connection #' + d.cid + d.title;
}

// The content for the popover on the connections
function popoverConnectionContent(d,i) {
  var content = '';
  var array = [{label: 'Type', value: d.type}, {label: 'Source Chromosome', value: d.sourceChromosome}, {label: 'Source Interval', value: Math.abs(d.source) + (d.source > 0 ? ' (tail)' : ' (head)')}, {label: 'Source Point', value: d3.format(',')(d.sourcePoint)},
   {label: 'Source Jabba', value: (d.sourceJabba)}, {label: 'Sink Chromosome', value: d.sinkChromosome}, {label: 'Sink Interval', value: Math.abs(d.sink) + (d.sink > 0 ? ' (tail)' : ' (head)')},
   {label: 'Sink Point', value: d3.format(',')(d.sinkPoint)}, {label: 'Sink Jabba', value: (d.sinkJabba)}, {label: 'Distance', value: (d.distance ? d3.format(',')(d.distance) : '-')}];
  if (!d.source) {
    array = [{label: 'Sink Chromosome', value: d.sinkChromosome}, {label: 'Sink Interval', value: Math.abs(d.sink) + (d.sink > 0 ? ' (tail)' : ' (head)')},
    {label: 'Sink Point', value: d3.format(',')(d.sinkPoint)}, {label: 'Sink Jabba', value: (d.sinkJabba)}, {label: 'Distance', value: (d.distance ? d3.format(',')(d.distance) : '-')}];
  }
  if (!d.sink) {
    array = [{label: 'Type', value: d.type}, {label: 'Source Chromosome', value: d.sourceChromosome}, {label: 'Source Interval', value: Math.abs(d.source) + (d.source > 0 ? ' (tail)' : ' (head)')}, {label: 'Source Point', value: d3.format(',')(d.sourcePoint)},
    {label: 'Source Jabba', value: (d.sourceJabba)}, {label: 'Distance', value: (d.distance ? d3.format(',')(d.distance) : '-')}];
  }
  array.forEach(function(e,j) {
     content += '<tr><td class="table-label" align="left" width="200" valign="top"><strong>' + e.label + ':</strong></td><td class="table-value" width="100" align="right" valign="top">' + e.value + '</td></tr>';
   });
  return '<div class="row"><div class="col-lg-12"><table width="0" border="0" align="left" cellpadding="0" cellspacing="0"><tbody>' + content + '</tbody></table></div></div>';
}

// The array of points forming the connections between its endpoints
function calculateConnectorEndpoints(yScale, record, connector, chromosome) {
  var points = [];
  record.sourceJabba = connector.source.y;
  record.sinkJabba = connector.sink.y;
  record.sourceChromosome = connector.source.chromosome;
  record.sinkChromosome = connector.sink.chromosome;
  record.sourcePoint = (connector.connection.source < 0) ?  connector.source.startPoint : connector.source.endPoint;
  record.sinkPoint   = (connector.connection.sink < 0)   ?  connector.sink.startPoint   : connector.sink.endPoint;
  record.distance = Math.abs(record.sinkPoint - record.sourcePoint);
  var origin = d3.min([record.sourcePoint, record.sinkPoint]);
  var target = d3.max([record.sourcePoint, record.sinkPoint]);
  var originSign = (origin === record.sourcePoint) ? connector.connection.source : connector.connection.sink;
  var targetSign = (target === record.sourcePoint) ? connector.connection.source : connector.connection.sink;
  var originY = (origin === record.sourcePoint) ? Math.abs(connector.source.y) : Math.abs(connector.sink.y);
  var targetY = (target === record.sourcePoint) ? Math.abs(connector.source.y) : Math.abs(connector.sink.y);
  var midPointX = 0.5 * chromosome.scale(origin) + 0.5 * chromosome.scale(target);
  var midPointY = 0.5 * originY + 0.5 * targetY;
  if (record.type === 'ALT' ) {
    if (Math.abs(connector.source.y) === Math.abs(connector.sink.y)) {
      points = [
              [chromosome.scale(origin), yScale(originY)],
              [d3.min([chromosome.scale(origin) + Math.sign(originSign) * 5, 0.95 * midPointX]), yScale(originY)],
              [d3.min([chromosome.scale(origin) + Math.sign(originSign) * 25, 0.95 * midPointX]), yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
              [midPointX, yScale((midPointY + (midPointY < 10 ? 0.75 : 10 )))],
              [d3.max([chromosome.scale(target) + Math.sign(targetSign) * 25, 1.05 * midPointX]), yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
              [d3.max([chromosome.scale(target) + Math.sign(targetSign) * 5, 1.05 * midPointX]), yScale(targetY)],
              [chromosome.scale(target), yScale(targetY)]];
    } else {
      points = [
              [chromosome.scale(origin), yScale(originY)],
              [chromosome.scale(origin) + Math.sign(originSign) * 5, yScale(originY)],
              [chromosome.scale(origin) + Math.sign(originSign) * 25, yScale((originY + Math.sign(targetY - originY) * (originY < 10 ? 0.25 : 5 )))],
              [chromosome.scale(target) + Math.sign(targetSign) * 25, yScale((targetY - Math.sign(targetY - originY) * (targetY < 10 ? 0.25 : 5 )))],
              [chromosome.scale(target) + Math.sign(targetSign) * 5, yScale(targetY)],
              [chromosome.scale(target), yScale(targetY)]];
    }
  } else {
    points = [
            [chromosome.scale(origin), yScale(originY)],
            [chromosome.scale(target), yScale(targetY)]];
  }
  return points;
}

// The array of points forming the connections between its endpoints in different chromosomes
function calculateInterConnectorEndpoints(yScale, record, connector, panelsArray) {
  var points = [];
  var sourceChromosome = panelsArray.filter(function(d,i) { return d.chromosome === connector.source.chromosome })[0];
  var sinkChromosome = panelsArray.filter(function(d,i) { return d.chromosome === connector.sink.chromosome })[0];

  record.sourceJabba = connector.source.y;
  record.sinkJabba = connector.sink.y;
  record.sourceChromosome = sourceChromosome.chromosome;
  record.sinkChromosome = sinkChromosome.chromosome;
  record.sourcePoint = (connector.connection.source < 0) ?  connector.source.startPoint : connector.source.endPoint;
  record.sinkPoint   = (connector.connection.sink < 0)   ?  connector.sink.startPoint   : connector.sink.endPoint;

  var origin = d3.min([sourceChromosome.panelScale(record.sourcePoint), sinkChromosome.panelScale(record.sinkPoint)]);
  var target = d3.max([sourceChromosome.panelScale(record.sourcePoint), sinkChromosome.panelScale(record.sinkPoint)]);
  var originSign = (origin === sourceChromosome.panelScale(record.sourcePoint)) ? connector.connection.source : connector.connection.sink;
  var targetSign = (target === sourceChromosome.panelScale(record.sourcePoint)) ? connector.connection.source : connector.connection.sink;
  var originY = (origin === sourceChromosome.panelScale(record.sourcePoint)) ? Math.abs(connector.source.y) : Math.abs(connector.sink.y);
  var targetY = (target === sourceChromosome.panelScale(record.sourcePoint)) ? Math.abs(connector.source.y) : Math.abs(connector.sink.y);
  var midPointX = 0.5 * origin + 0.5 * target;
  var midPointY = 0.5 * originY + 0.5 * targetY;

  if (record.type === 'ALT' ) {
    if (Math.abs(connector.source.y) === Math.abs(connector.sink.y)) {
      points = [
              [origin, yScale(originY)],
              [d3.min([origin + Math.sign(originSign) * 5, 0.95 * midPointX]), yScale(originY)],
              [d3.min([origin + Math.sign(originSign) * 25, 0.95 * midPointX]), yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
              [midPointX, yScale((midPointY + (midPointY < 10 ? 0.75 : 10 )))],
              [d3.max([target + Math.sign(targetSign) * 25, 1.05 * midPointX]), yScale((midPointY + (midPointY < 10 ? 0.5 : 5 )))],
              [d3.max([target + Math.sign(targetSign) * 5, 1.05 * midPointX]), yScale(targetY)],
              [target, yScale(targetY)]];
    } else {
      points = [
              [origin, yScale(originY)],
              [origin + Math.sign(originSign) * 5, yScale(originY)],
              [origin + Math.sign(originSign) * 25, yScale((originY + Math.sign(targetY - originY) * (originY < 10 ? 0.25 : 5 )))],
              [target + Math.sign(targetSign) * 25, yScale((targetY - Math.sign(targetY - originY) * (targetY < 10 ? 0.25 : 5 )))],
              [target + Math.sign(targetSign) * 5, yScale(targetY)],
              [target, yScale(targetY)]];
    }
  } else {
    points = [
            [origin, yScale(originY)],
            [target, yScale(targetY)]];
  }
  return points;
}

// The array of points forming the loose connections with one endpoint missing
function calculateLooseConnectorEndpoints(yScale, record, connector, chromosome) {
  if (!record.source && !record.sink) return;
  var touchpointX, touchpointY, touchpointSign;
  if (record.source && !record.sink) {
    record.sourceJabba = connector.source.y;
    record.sourceChromosome = connector.source.chromosome;
    record.sourcePoint = connector.connection.source > 0 ? connector.source.endPoint : connector.source.startPoint;
    touchpointX = record.sourcePoint;
    touchpointY = connector.source.y;
    touchpointSign = Math.sign(connector.connection.source);
  }
  if (!record.source && record.sink) {
    record.sinkJabba = connector.sink.y;
    record.sinkChromosome = connector.sink.chromosome;
    record.sinkPoint = connector.connection.sink > 0 ? connector.sink.endPoint : connector.sink.startPoint;
    touchpointX = record.sinkPoint;
    touchpointY = connector.sink.y;
    touchpointSign = Math.sign(connector.connection.sink);
  }
  return [
    [chromosome.scale(touchpointX), yScale(touchpointY)],
    [chromosome.scale(touchpointX) + touchpointSign * 15, yScale(touchpointY + (touchpointY < 10 ? 0.25 : 5 ))],
    [chromosome.scale(touchpointX) + touchpointSign * 5, yScale(touchpointY + (touchpointY < 10 ? 0.75 : 5 ))]];
}

// The array of points forming the connections with the other end in another chromosome
function calculateLocalInterConnectorEndpoints(yScale, record, connector, chromosomeObject) {
  var touchpointX, touchpointY, touchpointSign;
  if (connector.source.chromosome === chromosomeObject.chromosome) {
    record.sourceJabba = connector.source.y;
    record.sourceChromosome = connector.source.chromosome;
    record.sourcePoint = connector.connection.source > 0 ? connector.source.endPoint : connector.source.startPoint;
    touchpointX = record.sourcePoint;
    touchpointY = connector.source.y;
    touchpointSign = Math.sign(connector.connection.source);
  }
  if (connector.sink.chromosome === chromosomeObject.chromosome) {
    record.sinkJabba = connector.sink.y;
    record.sinkChromosome = connector.sink.chromosome;
    record.sinkPoint = connector.connection.sink > 0 ? connector.sink.endPoint : connector.sink.startPoint;
    touchpointX = record.sinkPoint;
    touchpointY = connector.sink.y;
    touchpointSign = Math.sign(connector.connection.sink);
  }
  return [
    [chromosomeObject.scale(touchpointX), yScale(touchpointY)],
    [chromosomeObject.scale(touchpointX) + touchpointSign * 15, yScale(touchpointY - (touchpointY < 10 ? 0.25 : 5 ))],
    [chromosomeObject.scale(touchpointX) - touchpointSign * 5, yScale(touchpointY - (touchpointY < 10 ? 0.75 : 5 ))]];
}
