// Add clickable logo to sidebar
$(document).ready(function() {
  var sidebar = $('#pkgdown-sidebar');
  if (sidebar.length) {
    var logoLink = $('<a>')
      .attr('href', 'https://alamarbio.com/')
      .attr('target', '_blank')
      .attr('rel', 'noopener noreferrer')
      .addClass('logo-link')
      .attr('title', 'Visit Alamar Biosciences');

    sidebar.prepend(logoLink);
  }
});
