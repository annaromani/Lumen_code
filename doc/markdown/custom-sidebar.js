// Custom Hierarchical Sidebar Navigation for FORD Documentation
// Generates collapsible menus from markdown heading structure

document.addEventListener('DOMContentLoaded', function() {
  const sidebarNav = document.querySelector('.nav.nav-pills.flex-column');
  const mainContent = document.getElementById('text');
  
  if (!sidebarNav || !mainContent) return;
  
  // Map of section titles to links (from index.md structure)
  const docStructure = {
    'THEORY & FUNDAMENTALS': [
      { title: 'Equilibrium Green\'s Functions', links: [
        { text: 'YAMBO Theory Documentation', href: 'YAMBO_THEORY_DOCUMENTATION.html' },
        { text: 'Theory Code Examples', href: 'YAMBO_THEORY_CODE_EXAMPLES.html' }
      ]},
      { title: 'Real-Time & Non-Equilibrium Systems', links: [
        { text: 'Real-Time TDDFT Theory', href: 'YAMBO_RT_THEORY_DOCUMENTATION.html' },
        { text: 'NEGF Theory Guide', href: 'YAMBO_RT_NEGF_THEORY.html' },
        { text: 'Screening in RT Systems', href: 'YAMBO_RT_SCREENING_THEORY.html' },
        { text: 'Kadanoff-Baym Theory & Implementation', href: 'KADANOFF_BAYM_IMPLEMENTATION.html' }
      ]},
      { title: 'Exciton-Phonon Coupling', links: [
        { text: 'Exciton-Phonon Theory', href: 'EXCPH_THEORY.html' },
        { text: 'Phonon-Assisted Transitions', href: 'EXCPH_IMPLEMENTATION.html' },
        { text: 'Exciton Lifetimes', href: 'EXCPH_IMPLEMENTATION.html' }
      ]},
      { title: 'Quick References', links: [
        { text: 'Formula Quick Reference', href: 'YAMBO_FORMULAS_QUICK_REFERENCE.html' },
        { text: 'NEGF Implementation Quick Reference', href: 'YAMBO_RT_NEGF_QUICK_REFERENCE.html' }
      ]}
    ],
    'IMPLEMENTATION & ARCHITECTURE': [
      { title: 'Exciton-Phonon Coupling', links: [
        { text: 'Exciton-Phonon Implementation', href: 'EXCPH_IMPLEMENTATION.html' },
        { text: 'Matrix Element Computation', href: 'EXCPH_IMPLEMENTATION.html' },
        { text: 'Optical Properties Calculation', href: 'EXCPH_IMPLEMENTATION.html' },
        { text: 'Lifetime Calculations', href: 'EXCPH_IMPLEMENTATION.html' }
      ]},
      { title: 'Real-Time Simulations', links: [
        { text: 'RT NEGF Implementation Analysis', href: 'YAMBO_RT_NEGF_IMPLEMENTATION.html' },
        { text: 'RT Screening Implementation', href: 'YAMBO_RT_SCREENING_IMPLEMENTATION.html' },
        { text: 'Two-Time Green\'s Function Reconstruction', href: 'TWO_TIME_GREEN_FUNCTION_RECONSTRUCTION.html' },
        { text: 'Two-Time Summary', href: 'TWO_TIME_RECONSTRUCTION_SUMMARY.html' }
      ]},
      { title: 'System Architecture', links: [
        { text: 'Code Flow & Architecture', href: 'YAMBO_CODE_FLOW_ARCHITECTURE.html' },
        { text: 'I/O System Architecture', href: 'IO_System_Architecture.html' },
        { text: 'I/O Quick Reference', href: 'IO_Quick_Reference.html' }
      ]},
      { title: 'Memory & Allocation', links: [
        { text: 'Memory Allocation System', href: 'YAMBO_ALLOC_DOCUMENTATION.html' },
        { text: 'Allocation Architecture', href: 'YAMBO_ALLOC_ARCHITECTURE.html' },
        { text: 'Allocation Examples', href: 'YAMBO_ALLOC_EXAMPLES.html' },
        { text: 'Allocation Quick Reference', href: 'YAMBO_ALLOC_QUICK_REFERENCE.html' }
      ]},
      { title: 'User Configuration', links: [
        { text: 'Input Parameters Reference', href: 'INPUT_FLAGS.html' }
      ]}
    ],
    'VISUAL GUIDES & FLOWCHARTS': [
      { title: 'Real-Time Dynamics', links: [
        { text: 'Real-Time Propagation Visual Guide', href: 'YAMBO_RT_VISUAL_GUIDE.html' },
        { text: 'YPP Real-Time Flowcharts', href: 'YPP_RT_FLOWCHARTS.html' },
        { text: 'YPP Real-Time Documentation', href: 'YPP_RT_DOCUMENTATION.html' }
      ]},
      { title: 'Theory Visualization', links: [
        { text: 'Theory Visual Guide', href: 'YAMBO_THEORY_VISUAL_GUIDE.html' }
      ]},
      { title: 'System Diagrams', links: [
        { text: 'I/O System Flowchart', href: 'IO_System_Flowchart.html' }
      ]}
    ]
  };
  
  // Clear existing sidebar content
  sidebarNav.innerHTML = '';
  
  // Build hierarchical structure
  let html = '';
  let sectionIndex = 0;
  
  for (const [mainSection, subsections] of Object.entries(docStructure)) {
    const collapseId = `collapse-section-${sectionIndex}`;
    const isFirst = sectionIndex === 0;
    
    html += `
      <div class="nav-section mb-2">
        <button class="btn btn-sm btn-outline-secondary w-100 text-start fw-bold" 
                type="button" data-bs-toggle="collapse" data-bs-target="#${collapseId}" 
                aria-expanded="${isFirst ? 'true' : 'false'}" aria-controls="${collapseId}">
          <i class="fas fa-chevron-down me-2"></i>${mainSection}
        </button>
        <div class="collapse ${isFirst ? 'show' : ''}" id="${collapseId}">
          <nav class="nav flex-column ms-2 mt-2">
    `;
    
    for (const subsection of subsections) {
      html += `
        <div class="nav-subsection mb-2">
          <small class="text-muted fw-bold d-block mb-1">${subsection.title}</small>
          <nav class="nav flex-column ms-2">
      `;
      
      for (const link of subsection.links) {
        html += `
          <a class="nav-link nav-link-doc pe-1" href="${link.href}" style="font-size: 0.9rem;">
            ${link.text}
          </a>
        `;
      }
      
      html += `
          </nav>
        </div>
      `;
    }
    
    html += `
          </nav>
        </div>
      </div>
    `;
    
    sectionIndex++;
  }
  
  // Add special link for main page
  html += `
    <hr class="my-2">
    <a class="nav-link fw-bold" href="mainpage.html" style="font-size: 0.95rem;">
      <i class="fas fa-home me-2"></i>Main Documentation
    </a>
  `;
  
  sidebarNav.innerHTML = html;
  
  // Highlight current page in sidebar
  const currentPage = window.location.pathname.split('/').pop();
  document.querySelectorAll('.nav-link-doc').forEach(link => {
    if (link.getAttribute('href') === currentPage) {
      link.classList.add('active', 'bg-light', 'text-dark');
      // Expand parent section
      const section = link.closest('.collapse');
      if (section) {
        const bsCollapse = new bootstrap.Collapse(section, { toggle: true });
      }
    }
  });
  
});
