import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
import glob
import base64
from PIL import Image
from io import BytesIO

# Configure page
st.set_page_config(
    page_title="RNA-Seq Analysis Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Custom CSS
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #4257b2;
        margin-bottom: 1rem;
    }
    .sub-header {
        font-size: 1.8rem;
        color: #3c9d9b;
        margin-bottom: 0.5rem;
        margin-top: 1rem;
    }
    .section-header {
        font-size: 1.5rem;
        color: #2c3e50;
        margin-top: 1rem;
        margin-bottom: 0.5rem;
    }
    .info-text {
        font-size: 1rem;
        color: #505050;
    }
    .highlight {
        background-color: #f0f2f6;
        padding: 1.5rem;
        border-radius: 0.5rem;
    }
    .stProgress .st-bo {
        background-color: #4257b2;
    }
</style>
""", unsafe_allow_html=True)

# Path to RNA-seq analysis results
# This would need to be configured by the user
@st.cache_data
def get_result_path():
    return st.session_state.get('result_path', None)

class RNASeqDashboard:
    def __init__(self):
        self.result_path = None
        self.samples = []
        self.contrasts = []
        self.has_loaded_data = False
    
    def load_project(self, project_path):
        """Load RNA-seq project data from the specified path."""
        self.result_path = project_path
        self.has_loaded_data = True
        
        # Load samples
        self.samples = []
        alignment_path = os.path.join(project_path, "03_alignment")
        if os.path.exists(alignment_path):
            sample_dirs = [d for d in os.listdir(alignment_path) 
                          if os.path.isdir(os.path.join(alignment_path, d))]
            self.samples = sample_dirs
        
        # Load contrasts
        self.contrasts = []
        de_path = os.path.join(project_path, "06_differential_expression")
        if os.path.exists(de_path):
            de_files = glob.glob(os.path.join(de_path, "*_results.csv"))
            self.contrasts = [os.path.basename(f).replace("_results.csv", "") for f in de_files]
        
        st.session_state['result_path'] = project_path
        st.session_state['samples'] = self.samples
        st.session_state['contrasts'] = self.contrasts
    
    def load_count_matrix(self):
        """Load the gene count matrix."""
        count_matrix_path = os.path.join(self.result_path, "04_counts", "gene_count_matrix.txt")
        if os.path.exists(count_matrix_path):
            try:
                return pd.read_csv(count_matrix_path, sep='\t', index_col=0)
            except Exception as e:
                st.error(f"Error loading count matrix: {str(e)}")
                return None
        else:
            return None
    
    def load_normalized_counts(self):
        """Load normalized expression data."""
        norm_counts_path = os.path.join(self.result_path, "05_expression", "normalized_counts.txt")
        if os.path.exists(norm_counts_path):
            try:
                return pd.read_csv(norm_counts_path, sep='\t', index_col=0)
            except Exception as e:
                st.error(f"Error loading normalized counts: {str(e)}")
                return None
        else:
            return None
    
    def load_de_results(self, contrast):
        """Load differential expression results for a specific contrast."""
        de_results_path = os.path.join(self.result_path, "06_differential_expression", f"{contrast}_results.csv")
        if os.path.exists(de_results_path):
            try:
                return pd.read_csv(de_results_path, index_col=0)
            except Exception as e:
                st.error(f"Error loading DE results: {str(e)}")
                return None
        else:
            return None
    
    def load_pca_results(self):
        """Load PCA results."""
        pca_path = os.path.join(self.result_path, "05_expression", "pca_results.txt")
        if os.path.exists(pca_path):
            try:
                return pd.read_csv(pca_path, sep='\t', index_col=0)
            except Exception as e:
                st.error(f"Error loading PCA results: {str(e)}")
                return None
        else:
            return None
    
    def render_dashboard(self):
        """Render the main dashboard interface."""
        # Header
        st.markdown('<div class="main-header">RNA-Seq Analysis Dashboard</div>', unsafe_allow_html=True)
        
        # Sidebar for project selection and navigation
        with st.sidebar:
            st.markdown("## Project Settings")
            
            # Project path input
            project_path = st.text_input("Enter the path to RNA-seq results:", value=get_result_path() or "")
            
            if st.button("Load Project"):
                if os.path.exists(project_path):
                    self.load_project(project_path)
                    st.success(f"Loaded project from {project_path}")
                else:
                    st.error("Invalid path. Please enter a valid directory path.")
            
            if self.has_loaded_data or get_result_path():
                # Navigation
                st.markdown("## Navigation")
                page = st.selectbox(
                    "Select a section to view:",
                    ["Project Overview", "Quality Control", "Expression Analysis", 
                     "Differential Expression", "Functional Enrichment"]
                )
                
                # Display sample and contrast info
                if 'samples' in st.session_state:
                    st.markdown("### Samples")
                    st.write(f"Total samples: {len(st.session_state['samples'])}")
                    with st.expander("View samples"):
                        for sample in st.session_state['samples']:
                            st.write(f"- {sample}")
                
                if 'contrasts' in st.session_state:
                    st.markdown("### Contrasts")
                    st.write(f"Total contrasts: {len(st.session_state['contrasts'])}")
                    with st.expander("View contrasts"):
                        for contrast in st.session_state['contrasts']:
                            st.write(f"- {contrast}")
            
            # About section
            with st.expander("About this dashboard"):
                st.markdown("""
                This dashboard visualizes results from the Advanced RNA-Seq Pipeline (ARSP).
                Navigate through different sections to explore your RNA-seq analysis results.
                
                **Features:**
                - Quality control visualization
                - Expression pattern analysis
                - Differential expression exploration
                - Functional enrichment visualization
                
                For more information, visit the GitHub repository.
                """)
        
        # Main content based on selected page
        if not self.has_loaded_data and not get_result_path():
            self._render_welcome_page()
        else:
            if get_result_path() and not self.has_loaded_data:
                self.load_project(get_result_path())
            
            if page == "Project Overview":
                self._render_project_overview()
            elif page == "Quality Control":
                self._render_quality_control()
            elif page == "Expression Analysis":
                self._render_expression_analysis()
            elif page == "Differential Expression":
                self._render_differential_expression()
            elif page == "Functional Enrichment":
                self._render_functional_enrichment()
    
    def _render_welcome_page(self):
        """Render the welcome page when no project is loaded."""
        st.markdown('<div class="highlight">', unsafe_allow_html=True)
        st.markdown("""
        # Welcome to the RNA-Seq Analysis Dashboard 🧬
        
        This interactive dashboard helps you explore and visualize your RNA-seq analysis results.
        
        ## To get started:
        1. Enter the path to your RNA-seq analysis results directory in the sidebar
        2. Click "Load Project" to load your data
        3. Navigate through different sections to explore your results
        
        ## Features
        - **Project Overview**: See summary statistics of your analysis
        - **Quality Control**: Explore read quality and alignment metrics
        - **Expression Analysis**: Analyze gene expression patterns and sample relationships
        - **Differential Expression**: Investigate differentially expressed genes
        - **Functional Enrichment**: Explore pathway and functional enrichment results
        
        ## Need demo data?
        If you don't have RNA-seq results yet, you can run the Advanced RNA-Seq Pipeline (ARSP) on your data first.
        """)
        st.markdown('</div>', unsafe_allow_html=True)
        
        # Example layout
        col1, col2 = st.columns(2)
        with col1:
            st.markdown('<div class="section-header">RNA-Seq Workflow</div>', unsafe_allow_html=True)
            st.image("https://raw.githubusercontent.com/yourusername/rnaseq-pipeline/main/docs/images/workflow.png", 
                    caption="RNA-Seq Analysis Workflow", use_column_width=True)
        
        with col2:
            st.markdown('<div class="section-header">Dashboard Features</div>', unsafe_allow_html=True)
            st.markdown("""
            - Interactive visualizations
            - Sample comparison tools
            - Gene expression exploration
            - Customizable plots
            - Export capabilities
            """)
    
    def _render_project_overview(self):
        """Render the project overview page."""
        st.markdown('<div class="sub-header">Project Overview</div>', unsafe_allow_html=True)
        
        # Load project summary
        summary_path = os.path.join(self.result_path, "08_reports", "pipeline_summary.txt")
        if os.path.exists(summary_path):
            with open(summary_path, 'r') as f:
                summary_text = f.read()
            
            # Display summary information
            st.markdown('<div class="highlight">', unsafe_allow_html=True)
            st.text(summary_text)
            st.markdown('</div>', unsafe_allow_html=True)
        
        # Project statistics
        st.markdown('<div class="section-header">Project Statistics</div>', unsafe_allow_html=True)
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Samples", len(self.samples))
        
        with col2:
            st.metric("Total Genes", self._get_gene_count())
        
        with col3:
            st.metric("DE Contrasts", len(self.contrasts))
        
        with col4:
            st.metric("Significant DE Genes", self._get_significant_de_gene_count())
        
        # Display sample information
        st.markdown('<div class="section-header">Sample Information</div>', unsafe_allow_html=True)
        
        sample_stats = self._get_sample_statistics()
        if sample_stats is not None:
            st.dataframe(sample_stats)
        
        # Display analysis completion status
        st.markdown('<div class="section-header">Analysis Completion</div>', unsafe_allow_html=True)
        
        completion_status = {
            "Quality Control": os.path.exists(os.path.join(self.result_path, "01_fastqc")),
            "Read Trimming": os.path.exists(os.path.join(self.result_path, "02_trimmed")),
            "Read Alignment": os.path.exists(os.path.join(self.result_path, "03_alignment")),
            "Expression Quantification": os.path.exists(os.path.join(self.result_path, "04_counts")),
            "Expression Analysis": os.path.exists(os.path.join(self.result_path, "05_expression")),
            "Differential Expression": os.path.exists(os.path.join(self.result_path, "06_differential_expression")),
            "Functional Enrichment": os.path.exists(os.path.join(self.result_path, "06_differential_expression", "functional")),
            "Report Generation": os.path.exists(os.path.join(self.result_path, "08_reports"))
        }
        
        for step, completed in completion_status.items():
            st.write(f"{step}: {'✅' if completed else '❌'}")
    
    def _render_quality_control(self):
        """Render the quality control visualization page."""
        st.markdown('<div class="sub-header">Quality Control</div>', unsafe_allow_html=True)
        
        # Display MultiQC results if available
        raw_multiqc = os.path.join(self.result_path, "01_fastqc", "multiqc_report.html")
        trimmed_multiqc = os.path.join(self.result_path, "02_trimmed", "fastqc", "multiqc_report.html")
        alignment_multiqc = os.path.join(self.result_path, "03_alignment", "multiqc_report.html")
        
        # Options for what to view
        qc_option = st.radio(
            "Select QC results to view:",
            ["Raw Reads QC", "Trimmed Reads QC", "Alignment QC", "Read Mapping Statistics"]
        )
        
        if qc_option == "Raw Reads QC":
            if os.path.exists(raw_multiqc):
                st.markdown('<div class="section-header">Raw Reads Quality Control</div>', unsafe_allow_html=True)
                self._embed_html_file(raw_multiqc)
            else:
                self._display_qc_plots("01_fastqc")
        
        elif qc_option == "Trimmed Reads QC":
            if os.path.exists(trimmed_multiqc):
                st.markdown('<div class="section-header">Trimmed Reads Quality Control</div>', unsafe_allow_html=True)
                self._embed_html_file(trimmed_multiqc)
            else:
                self._display_qc_plots("02_trimmed/fastqc")
        
        elif qc_option == "Alignment QC":
            if os.path.exists(alignment_multiqc):
                st.markdown('<div class="section-header">Alignment Quality Control</div>', unsafe_allow_html=True)
                self._embed_html_file(alignment_multiqc)
            else:
                st.markdown('<div class="info-text">MultiQC alignment report not found.</div>', unsafe_allow_html=True)
        
        elif qc_option == "Read Mapping Statistics":
            self._display_alignment_stats()
    
    def _display_alignment_stats(self):
        """Display alignment statistics for all samples."""
        st.markdown('<div class="section-header">Read Mapping Statistics</div>', unsafe_allow_html=True)
        
        # Collect alignment statistics
        alignment_stats = {}
        for sample in self.samples:
            flagstat_file = os.path.join(self.result_path, "03_alignment", sample, f"{sample}_flagstat.txt")
            if os.path.exists(flagstat_file):
                with open(flagstat_file, 'r') as f:
                    stats = f.readlines()
                
                # Parse flagstat output
                total_reads = int(stats[0].split()[0])
                mapped_reads = int(stats[4].split()[0])
                mapping_rate = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0
                
                alignment_stats[sample] = {
                    'Total Reads': total_reads,
                    'Mapped Reads': mapped_reads,
                    'Mapping Rate (%)': round(mapping_rate, 2)
                }
        
        if alignment_stats:
            # Convert to DataFrame
            stats_df = pd.DataFrame.from_dict(alignment_stats, orient='index')
            
            # Display table
            st.dataframe(stats_df)
            
            # Display bar chart of mapping rates
            fig = px.bar(
                stats_df, 
                y='Mapping Rate (%)', 
                title='Read Mapping Rates by Sample',
                height=500
            )
            fig.update_layout(
                xaxis_title='Sample',
                yaxis_title='Mapping Rate (%)',
                yaxis_range=[0, 100]
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.markdown('<div class="info-text">No alignment statistics found.</div>', unsafe_allow_html=True)
    
    def _embed_html_file(self, file_path):
        """Embed an HTML file in the dashboard."""
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html_content = f.read()
            
            # Use an iframe to display the HTML content
            st.components.v1.html(html_content, height=600, scrolling=True)
            
            # Provide download link
            st.download_button(
                label="Download HTML Report",
                data=html_content,
                file_name=os.path.basename(file_path),
                mime="text/html"
            )
        except Exception as e:
            st.error(f"Error loading HTML file: {str(e)}")
    
    def _display_qc_plots(self, qc_dir_path):
        """Display QC plots from the specified directory."""
        qc_dir = os.path.join(self.result_path, qc_dir_path)
        
        if not os.path.exists(qc_dir):
            st.markdown('<div class="info-text">QC directory not found.</div>', unsafe_allow_html=True)
            return
        
        # Find all PNG files in the QC directory
        png_files = glob.glob(os.path.join(qc_dir, "**", "*.png"), recursive=True)
        
        if not png_files:
            st.markdown('<div class="info-text">No QC plots found.</div>', unsafe_allow_html=True)
            return
        
        # Group PNG files by sample
        png_by_sample = {}
        for png_file in png_files:
            filename = os.path.basename(png_file)
            sample = filename.split("_fastqc")[0]
            if sample not in png_by_sample:
                png_by_sample[sample] = []
            png_by_sample[sample].append(png_file)
        
        # Select sample to view
        selected_sample = st.selectbox("Select sample:", list(png_by_sample.keys()))
        
        if selected_sample:
            # Display plots for selected sample
            for png_file in png_by_sample[selected_sample]:
                try:
                    img = Image.open(png_file)
                    st.image(img, caption=os.path.basename(png_file), use_column_width=True)
                except Exception as e:
                    st.error(f"Error loading image {os.path.basename(png_file)}: {str(e)}")
    
    def _render_expression_analysis(self):
        """Render the expression analysis visualization page."""
        st.markdown('<div class="sub-header">Expression Analysis</div>', unsafe_allow_html=True)
        
        # Tabs for different expression analysis visualizations
        tabs = st.tabs(["Expression Overview", "PCA Analysis", "Sample Correlation", "Gene Expression"])
        
        with tabs[0]:  # Expression Overview
            self._render_expression_overview()
        
        with tabs[1]:  # PCA Analysis
            self._render_pca_analysis()
        
        with tabs[2]:  # Sample Correlation
            self._render_sample_correlation()
        
        with tabs[3]:  # Gene Expression
            self._render_gene_expression_explorer()
    
    def _render_expression_overview(self):
        """Render expression overview visualizations."""
        st.markdown('<div class="section-header">Expression Level Distribution</div>', unsafe_allow_html=True)
        
        # Load normalized counts
        norm_counts = self.load_normalized_counts()
        if norm_counts is None:
            st.markdown('<div class="info-text">Normalized expression data not found.</div>', unsafe_allow_html=True)
            return
        
        # Expression level distribution
        fig = px.box(
            pd.DataFrame(norm_counts.transpose()), 
            title='Expression Level Distribution by Sample', 
            height=600,
            log_y=True
        )
        fig.update_layout(
            xaxis_title='Sample',
            yaxis_title='Expression Level (log scale)',
            showlegend=False
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Expression statistics
        st.markdown('<div class="section-header">Expression Statistics</div>', unsafe_allow_html=True)
        
        # Calculate basic statistics
        expr_stats = pd.DataFrame({
            'Mean': norm_counts.mean(axis=1),
            'Median': norm_counts.median(axis=1),
            'Min': norm_counts.min(axis=1),
            'Max': norm_counts.max(axis=1),
            'Std Dev': norm_counts.std(axis=1)
        })
        
        # Display statistics for top expressed genes
        top_genes = expr_stats.sort_values('Mean', ascending=False).head(20)
        
        st.markdown("**Top 20 Highly Expressed Genes:**")
        st.dataframe(top_genes)
        
        # Plot distribution of expression values
        fig = px.histogram(
            expr_stats, 
            x='Mean', 
            title='Distribution of Mean Expression Values',
            nbins=50,
            log_x=True
        )
        fig.update_layout(
            xaxis_title='Mean Expression (log scale)',
            yaxis_title='Count'
        )
        st.plotly_chart(fig, use_container_width=True)
    
    def _render_pca_analysis(self):
        """Render PCA analysis visualizations."""
        st.markdown('<div class="section-header">Principal Component Analysis</div>', unsafe_allow_html=True)
        
        # Load PCA results
        pca_results = self.load_pca_results()
        if pca_results is None:
            # If PCA results don't exist, compute PCA from normalized counts
            norm_counts = self.load_normalized_counts()
            if norm_counts is None:
                st.markdown('<div class="info-text">Cannot perform PCA: No expression data found.</div>', unsafe_allow_html=True)
                return
            
            # Perform PCA
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import StandardScaler
            
            # Transpose so that samples are rows
            data = norm_counts.transpose()
            
            # Filter low-variance genes
            data = data.loc[:, data.var() > 0.1]
            
            # Standardize data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(data)
            
            # Apply PCA
            n_components = min(3, scaled_data.shape[0], scaled_data.shape[1])
            pca = PCA(n_components=n_components)
            pca_data = pca.fit_transform(scaled_data)
            
            # Create DataFrame with PCA results
            pca_results = pd.DataFrame(
                pca_data,
                index=data.index,
                columns=[f'PC{i+1}' for i in range(n_components)]
            )
            
            # Add variance explained
            variance_explained = pca.explained_variance_ratio_ * 100
            
            st.info("PCA results were computed from normalized counts.")
        else:
            # Load variance explained
            variance_path = os.path.join(self.result_path, "05_expression", "pca_variance.txt")
            if os.path.exists(variance_path):
                variance_df = pd.read_csv(variance_path, sep='\t', index_col=0)
                variance_explained = variance_df['Variance_explained'].values * 100
            else:
                variance_explained = None
        
        # Load sample metadata if available
        metadata = self._load_sample_metadata()
        
        # Create PCA plot
        if metadata is not None and set(pca_results.index).issubset(set(metadata.index)):
            # Merge PCA results with metadata
            plot_data = pca_results.copy()
            for col in metadata.columns:
                plot_data[col] = metadata.loc[plot_data.index, col].values
            
            # Create scatter plot with color by condition
            color_col = 'condition' if 'condition' in metadata.columns else None
            symbol_col = 'batch' if 'batch' in metadata.columns else None
            
            fig = px.scatter(
                plot_data,
                x='PC1',
                y='PC2',
                color=color_col,
                symbol=symbol_col,
                text=plot_data.index,
                title='PCA Plot',
                height=600
            )
        else:
            # Basic PCA plot without metadata
            fig = px.scatter(
                pca_results,
                x='PC1',
                y='PC2',
                text=pca_results.index,
                title='PCA Plot',
                height=600
            )
        
        # Add variance explained to axis labels
        if variance_explained is not None and len(variance_explained) >= 2:
            fig.update_layout(
                xaxis_title=f'PC1 ({variance_explained[0]:.1f}% variance)',
                yaxis_title=f'PC2 ({variance_explained[1]:.1f}% variance)'
            )
        
        # Show sample names on hover
        fig.update_traces(
            hovertemplate='<b>%{text}</b><br>PC1: %{x:.4f}<br>PC2: %{y:.4f}'
        )
        
        st.plotly_chart(fig, use_container_width=True)
        
        # If we have at least 3 PCs, add 3D plot option
        if 'PC3' in pca_results.columns:
            show_3d = st.checkbox("Show 3D PCA Plot")
            
            if show_3d:
                if metadata is not None and set(pca_results.index).issubset(set(metadata.index)):
                    # Merge PCA results with metadata
                    plot_data = pca_results.copy()
                    for col in metadata.columns:
                        plot_data[col] = metadata.loc[plot_data.index, col].values
                    
                    # Create 3D scatter plot with color by condition
                    color_col = 'condition' if 'condition' in metadata.columns else None
                    
                    fig_3d = px.scatter_3d(
                        plot_data,
                        x='PC1',
                        y='PC2',
                        z='PC3',
                        color=color_col,
                        text=plot_data.index,
                        title='3D PCA Plot',
                        height=700
                    )
                else:
                    # Basic 3D PCA plot without metadata
                    fig_3d = px.scatter_3d(
                        pca_results,
                        x='PC1',
                        y='PC2',
                        z='PC3',
                        text=pca_results.index,
                        title='3D PCA Plot',
                        height=700
                    )
                
                # Add variance explained to axis labels
                if variance_explained is not None and len(variance_explained) >= 3:
                    fig_3d.update_layout(
                        scene=dict(
                            xaxis_title=f'PC1 ({variance_explained[0]:.1f}% variance)',
                            yaxis_title=f'PC2 ({variance_explained[1]:.1f}% variance)',
                            zaxis_title=f'PC3 ({variance_explained[2]:.1f}% variance)'
                        )
                    )
                
                # Show sample names on hover
                fig_3d.update_traces(
                    hovertemplate='<b>%{text}</b><br>PC1: %{x:.4f}<br>PC2: %{y:.4f}<br>PC3: %{z:.4f}'
                )
                
                st.plotly_chart(fig_3d, use_container_width=True)
        
        # Variance explained
        if variance_explained is not None:
            st.markdown('<div class="section-header">Variance Explained by Principal Components</div>', unsafe_allow_html=True)
            
            # Create bar chart of variance explained
            fig = px.bar(
                x=[f'PC{i+1}' for i in range(len(variance_explained))],
                y=variance_explained,
                title='Variance Explained by Each Principal Component',
                height=400
            )
            fig.update_layout(
                xaxis_title='Principal Component',
                yaxis_title='Variance Explained (%)',
                yaxis_range=[0, max(100, max(variance_explained) * 1.1)]
            )
            st.plotly_chart(fig, use_container_width=True)
    
    def _render_sample_correlation(self):
        """Render sample correlation visualizations."""
        st.markdown('<div class="section-header">Sample Correlation Analysis</div>', unsafe_allow_html=True)
        
        # Load normalized counts
        norm_counts = self.load_normalized_counts()
        if norm_counts is None:
            st.markdown('<div class="info-text">Normalized expression data not found.</div>', unsafe_allow_html=True)
            return
        
        # Calculate correlation matrix
        corr_matrix = norm_counts.corr()
        
        # Create heatmap
        fig = px.imshow(
            corr_matrix,
            title='Sample Correlation Heatmap',
            color_continuous_scale='Viridis',
            zmin=0.5,
            zmax=1
        )
        fig.update_layout(
            height=600,
            xaxis_title='Sample',
            yaxis_title='Sample'
        )
        st.plotly_chart(fig, use_container_width=True)
        
        # Hierarchical clustering of samples
        st.markdown('<div class="section-header">Sample Clustering</div>', unsafe_allow_html=True)
        
        # Convert correlation to distance
        from scipy.spatial.distance import squareform
        from scipy.cluster.hierarchy import linkage, dendrogram
        
        # Calculate distance matrix (1 - correlation)
        distance_matrix = 1 - corr_matrix
        
        # Calculate linkage
        Z = linkage(squareform(distance_matrix), method='average')
        
        # Create dendrogram
        fig, ax = plt.subplots(figsize=(12, 6))
        dendrogram(
            Z,
            labels=corr_matrix.index,
            leaf_rotation=90,
            ax=ax
        )
        ax.set_title('Hierarchical Clustering of Samples')
        st.pyplot(fig)
    
    def _render_gene_expression_explorer(self):
        """Render gene expression explorer."""
        st.markdown('<div class="section-header">Gene Expression Explorer</div>', unsafe_allow_html=True)
        
        # Load normalized counts
        norm_counts = self.load_normalized_counts()
        if norm_counts is None:
            st.markdown('<div class="info-text">Normalized expression data not found.</div>', unsafe_allow_html=True)
            return
        
        # Gene selection
        gene_search = st.text_input("Search for gene (gene ID or partial name):")
        
        if gene_search:
            # Find matching genes
            matching_genes = [gene for gene in norm_counts.index if gene_search.upper() in gene.upper()]
            
            if matching_genes:
                if len(matching_genes) > 100:
                    st.warning(f"Found {len(matching_genes)} matching genes. Showing first 100.")
                    matching_genes = matching_genes[:100]
                
                selected_gene = st.selectbox("Select gene:", matching_genes)
                
                if selected_gene:
                    self._display_gene_expression(selected_gene, norm_counts)
            else:
                st.info(f"No genes found matching '{gene_search}'.")
        else:
            # Allow selection from top expressed genes
            top_genes = norm_counts.mean(axis=1).sort_values(ascending=False).head(100).index.tolist()
            st.markdown('<div class="info-text">Select from top 100 expressed genes:</div>', unsafe_allow_html=True)
            
            selected_gene = st.selectbox("Select gene:", top_genes)
            
            if selected_gene:
                self._display_gene_expression(selected_gene, norm_counts)
    
    def _display_gene_expression(self, gene, norm_counts):
        """Display expression data for a selected gene."""
        # Get expression values for the gene
        if gene not in norm_counts.index:
            st.error(f"Gene {gene} not found in expression data.")
            return
        
        expr_values = norm_counts.loc[gene]
        
        # Load sample metadata if available
        metadata = self._load_sample_metadata()
        
        # Create plot data
        if metadata is not None and set(expr_values.index).issubset(set(metadata.index)):
            # Merge expression with metadata
            plot_data = pd.DataFrame({'Expression': expr_values})
            for col in metadata.columns:
                plot_data[col] = metadata.loc[plot_data.index, col].values
            
            # Check if we have condition information
            if 'condition' in metadata.columns:
                # Group by condition
                grouped_data = plot_data.groupby('condition')['Expression'].agg(['mean', 'std']).reset_index()
                
                # Create bar chart with error bars
                fig = make_subplots(rows=1, cols=2, subplot_titles=["Expression by Condition", "Individual Samples"])
                
                # Add bar chart for condition means
                fig.add_trace(
                    go.Bar(
                        x=grouped_data['condition'],
                        y=grouped_data['mean'],
                        error_y=dict(type='data', array=grouped_data['std']),
                        name='Mean Expression'
                    ),
                    row=1, col=1
                )
                
                # Add scatter plot for individual samples
                fig.add_trace(
                    go.Scatter(
                        x=plot_data['condition'],
                        y=plot_data['Expression'],
                        mode='markers',
                        marker=dict(size=10),
                        name='Individual Samples'
                    ),
                    row=1, col=2
                )
                
                fig.update_layout(
                    title=f'Expression of {gene}',
                    height=500,
                    showlegend=False
                )
                
                st.plotly_chart(fig, use_container_width=True)
                
                # Statistical test if we have multiple conditions
                if len(grouped_data) > 1:
                    from scipy.stats import f_oneway
                    
                    # One-way ANOVA
                    groups = [plot_data[plot_data['condition'] == c]['Expression'].values 
                             for c in grouped_data['condition']]
                    f_stat, p_value = f_oneway(*groups)
                    
                    st.markdown(f"**One-way ANOVA:** F-statistic = {f_stat:.4f}, p-value = {p_value:.4e}")
                    
                    if p_value < 0.05:
                        st.markdown("**Result:** There is a significant difference in expression between conditions.")
                    else:
                        st.markdown("**Result:** No significant difference in expression between conditions.")
            else:
                # Simple bar chart without grouping
                fig = px.bar(
                    x=expr_values.index,
                    y=expr_values.values,
                    title=f'Expression of {gene}',
                    height=500
                )
                fig.update_layout(
                    xaxis_title='Sample',
                    yaxis_title='Expression'
                )
                st.plotly_chart(fig, use_container_width=True)
        else:
            # Simple bar chart without metadata
            fig = px.bar(
                x=expr_values.index,
                y=expr_values.values,
                title=f'Expression of {gene}',
                height=500
            )
            fig.update_layout(
                xaxis_title='Sample',
                yaxis_title='Expression'
            )
            st.plotly_chart(fig, use_container_width=True)
        
        # Display gene expression values
        st.markdown('<div class="section-header">Expression Values</div>', unsafe_allow_html=True)
        expr_df = pd.DataFrame({'Expression': expr_values})
        st.dataframe(expr_df)
        
        # Check if gene is DE in any contrast
        st.markdown('<div class="section-header">Differential Expression Status</div>', unsafe_allow_html=True)
        
        de_status = {}
        for contrast in self.contrasts:
            de_results = self.load_de_results(contrast)
            if de_results is not None and gene in de_results.index:
                # Get log2FC and adjusted p-value
                log2fc = de_results.loc[gene, 'log2FoldChange'] if 'log2FoldChange' in de_results.columns else de_results.loc[gene, 'logFC']
                padj = de_results.loc[gene, 'padj'] if 'padj' in de_results.columns else de_results.loc[gene, 'FDR']
                
                de_status[contrast] = {
                    'Log2FC': log2fc,
                    'Adj. P-value': padj,
                    'Significant': padj < 0.05 and abs(log2fc) > 1
                }
        
        if de_status:
            de_df = pd.DataFrame.from_dict(de_status, orient='index')
            
            # Add colored background for significant DE
            st.dataframe(de_df.style.apply(
                lambda row: ['background-color: #d4edda' if row['Significant'] else '' for _ in row],
                axis=1
            ))
        else:
            st.markdown('<div class="info-text">No differential expression results found for this gene.</div>', unsafe_allow_html=True)
    
    def _render_differential_expression(self):
        """Render differential expression visualizations."""
        st.markdown('<div class="sub-header">Differential Expression Analysis</div>', unsafe_allow_html=True)
        
        if not self.contrasts:
            st.markdown('<div class="info-text">No differential expression contrasts found.</div>', unsafe_allow_html=True)
            return
        
        # Select contrast to visualize
        selected_contrast = st.selectbox("Select contrast:", self.contrasts)
        
        if selected_contrast:
            # Load DE results
            de_results = self.load_de_results(selected_contrast)
            
            if de_results is None:
                st.markdown('<div class="info-text">Differential expression results not found for this contrast.</div>', unsafe_allow_html=True)
                return
            
            # Determine columns based on DE method (DESeq2 or edgeR)
            if 'log2FoldChange' in de_results.columns:
                # DESeq2 format
                fc_col = 'log2FoldChange'
                pval_col = 'pvalue'
                padj_col = 'padj'
            else:
                # edgeR format
                fc_col = 'logFC'
                pval_col = 'PValue'
                padj_col = 'FDR'
            
            # Add volcano plot
            st.markdown('<div class="section-header">Volcano Plot</div>', unsafe_allow_html=True)
            
            # Add significance and color
            de_results['significant'] = (de_results[padj_col] < 0.05) & (abs(de_results[fc_col]) > 1)
            de_results['regulation'] = 'Not Significant'
            de_results.loc[(de_results['significant']) & (de_results[fc_col] > 1), 'regulation'] = 'Up-regulated'
            de_results.loc[(de_results['significant']) & (de_results[fc_col] < -1), 'regulation'] = 'Down-regulated'
            
            # Create volcano plot
            fig = px.scatter(
                de_results,
                x=fc_col,
                y=-np.log10(de_results[pval_col]),
                color='regulation',
                color_discrete_map={
                    'Up-regulated': 'red',
                    'Down-regulated': 'blue',
                    'Not Significant': 'gray'
                },
                opacity=0.7,
                hover_name=de_results.index,
                hover_data={
                    fc_col: ':.3f',
                    pval_col: ':.2e',
                    padj_col: ':.2e'
                },
                title=f'Volcano Plot: {selected_contrast}'
            )
            
            # Add threshold lines
            fig.add_shape(
                type="line", line=dict(dash="dash"),
                x0=-1, y0=0, x1=-1, y1=int(-np.log10(de_results[pval_col].min())) + 1,
                line_color="black"
            )
            fig.add_shape(
                type="line", line=dict(dash="dash"),
                x0=1, y0=0, x1=1, y1=int(-np.log10(de_results[pval_col].min())) + 1,
                line_color="black"
            )
            fig.add_shape(
                type="line", line=dict(dash="dash"),
                x0=-max(abs(de_results[fc_col])) - 0.5, y0=-np.log10(0.05), 
                x1=max(abs(de_results[fc_col])) + 0.5, y1=-np.log10(0.05),
                line_color="black"
            )
            
            fig.update_layout(
                xaxis_title=f"{fc_col}",
                yaxis_title="-log10(p-value)",
                height=600
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Show MA plot if available
            ma_plot_path = os.path.join(self.result_path, "06_differential_expression", f"{selected_contrast}_MA_plot.png")
            if os.path.exists(ma_plot_path):
                st.markdown('<div class="section-header">MA Plot</div>', unsafe_allow_html=True)
                try:
                    img = Image.open(ma_plot_path)
                    st.image(img, use_column_width=True)
                except Exception as e:
                    st.error(f"Error loading MA plot: {str(e)}")
            
            # DE genes summary
            st.markdown('<div class="section-header">Differential Expression Summary</div>', unsafe_allow_html=True)
            
            # Count significant genes
            sig_up = ((de_results[padj_col] < 0.05) & (de_results[fc_col] > 1)).sum()
            sig_down = ((de_results[padj_col] < 0.05) & (de_results[fc_col] < -1)).sum()
            total_genes = len(de_results)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Total Genes", total_genes)
            col2.metric("Upregulated Genes", sig_up)
            col3.metric("Downregulated Genes", sig_down)
            
            # Top DE genes table
            st.markdown('<div class="section-header">Top Differentially Expressed Genes</div>', unsafe_allow_html=True)
            
            # Sort by adjusted p-value
            top_de = de_results.sort_values(padj_col).head(20)
            
            # Format for display
            display_df = top_de[[fc_col, pval_col, padj_col, 'regulation']].copy()
            display_df.columns = ['Log2 Fold Change', 'P-value', 'Adj. P-value', 'Regulation']
            
            # Format numeric columns
            display_df['P-value'] = display_df['P-value'].apply(lambda x: f"{x:.2e}")
            display_df['Adj. P-value'] = display_df['Adj. P-value'].apply(lambda x: f"{x:.2e}")
            display_df['Log2 Fold Change'] = display_df['Log2 Fold Change'].apply(lambda x: f"{x:.3f}")
            
            # Display table
            st.dataframe(display_df)
            
            # Interactive gene search
            st.markdown('<div class="section-header">Search for Specific Genes</div>', unsafe_allow_html=True)
            
            gene_search = st.text_input("Enter gene ID or partial name:")
            
            if gene_search:
                # Find matching genes
                matching_genes = de_results[de_results.index.str.contains(gene_search, case=False)]
                
                if not matching_genes.empty:
                    st.write(f"Found {len(matching_genes)} matching genes.")
                    st.dataframe(matching_genes[[fc_col, pval_col, padj_col, 'regulation']])
                else:
                    st.info(f"No genes found matching '{gene_search}'.")
            
            # Download results option
            st.markdown('<div class="section-header">Download Results</div>', unsafe_allow_html=True)
            
            # Create a CSV for download
            csv = de_results.to_csv()
            st.download_button(
                label="Download Complete DE Results",
                data=csv,
                file_name=f"{selected_contrast}_de_results.csv",
                mime="text/csv"
            )
    
    def _render_functional_enrichment(self):
        """Render functional enrichment visualizations."""
        st.markdown('<div class="sub-header">Functional Enrichment Analysis</div>', unsafe_allow_html=True)
        
        if not self.contrasts:
            st.markdown('<div class="info-text">No contrasts found for functional enrichment analysis.</div>', unsafe_allow_html=True)
            return
        
        # Select contrast to visualize
        selected_contrast = st.selectbox("Select contrast:", self.contrasts)
        
        if selected_contrast:
            # Path to functional enrichment results
            func_dir = os.path.join(self.result_path, "06_differential_expression", "functional", selected_contrast)
            
            if not os.path.exists(func_dir):
                st.markdown('<div class="info-text">Functional enrichment results not found for this contrast.</div>', unsafe_allow_html=True)
                return
            
            # Tabs for different types of enrichment
            tabs = st.tabs(["GO Biological Process", "GO Molecular Function", "GO Cellular Component", "KEGG Pathways"])
            
            with tabs[0]:  # GO BP
                self._display_enrichment_results(func_dir, "GO_BP")
            
            with tabs[1]:  # GO MF
                self._display_enrichment_results(func_dir, "GO_MF")
            
            with tabs[2]:  # GO CC
                self._display_enrichment_results(func_dir, "GO_CC")
            
            with tabs[3]:  # KEGG
                self._display_enrichment_results(func_dir, "KEGG")
    
    def _display_enrichment_results(self, func_dir, analysis_type):
        """Display enrichment results for a specific analysis type."""
        # Load enrichment results if available
        enrich_file = os.path.join(func_dir, f"{analysis_type}_enrichment.csv")
        
        if os.path.exists(enrich_file):
            try:
                enrich_df = pd.read_csv(enrich_file, index_col=0)
                
                # Display enrichment table
                st.markdown(f'<div class="section-header">{analysis_type.replace("_", " ")} Enrichment Results</div>', unsafe_allow_html=True)
                
                if not enrich_df.empty:
                    # Format p-values
                    if 'pvalue' in enrich_df.columns:
                        enrich_df['pvalue'] = enrich_df['pvalue'].apply(lambda x: f"{x:.2e}")
                    if 'p.adjust' in enrich_df.columns:
                        enrich_df['p.adjust'] = enrich_df['p.adjust'].apply(lambda x: f"{x:.2e}")
                    
                    st.dataframe(enrich_df.head(20))
                    
                    # Display visualization if available
                    barplot_file = os.path.join(func_dir, f"{analysis_type}_barplot.png")
                    dotplot_file = os.path.join(func_dir, f"{analysis_type}_dotplot.png")
                    
                    if os.path.exists(barplot_file):
                        st.markdown('<div class="section-header">Enrichment Barplot</div>', unsafe_allow_html=True)
                        try:
                            img = Image.open(barplot_file)
                            st.image(img, use_column_width=True)
                        except Exception as e:
                            st.error(f"Error loading barplot: {str(e)}")
                    
                    if os.path.exists(dotplot_file):
                        st.markdown('<div class="section-header">Enrichment Dotplot</div>', unsafe_allow_html=True)
                        try:
                            img = Image.open(dotplot_file)
                            st.image(img, use_column_width=True)
                        except Exception as e:
                            st.error(f"Error loading dotplot: {str(e)}")
                    
                    # KEGG pathway plots
                    if analysis_type == "KEGG" and 'ID' in enrich_df.columns:
                        pathway_files = glob.glob(os.path.join(func_dir, "KEGG_*.png"))
                        
                        if pathway_files:
                            st.markdown('<div class="section-header">KEGG Pathway Visualizations</div>', unsafe_allow_html=True)
                            
                            for pathway_file in pathway_files[:3]:  # Limit to first 3 pathways
                                pathway_name = os.path.basename(pathway_file).replace("KEGG_", "").replace(".png", "")
                                try:
                                    img = Image.open(pathway_file)
                                    st.image(img, caption=pathway_name, use_column_width=True)
                                except Exception as e:
                                    st.error(f"Error loading pathway image: {str(e)}")
                    
                    # Download option
                    csv = enrich_df.to_csv()
                    st.download_button(
                        label=f"Download {analysis_type} Enrichment Results",
                        data=csv,
                        file_name=f"{analysis_type}_enrichment.csv",
                        mime="text/csv"
                    )
                else:
                    st.info(f"No significant {analysis_type} enrichment results found.")
            except Exception as e:
                st.error(f"Error loading enrichment results: {str(e)}")
        else:
            st.info(f"No {analysis_type} enrichment results found.")
    
    def _get_gene_count(self):
        """Get the total number of genes in the count matrix."""
        count_matrix = self.load_count_matrix()
        if count_matrix is not None:
            return count_matrix.shape[0]
        return 0
    
    def _get_significant_de_gene_count(self):
        """Get the total number of significant DE genes across all contrasts."""
        total_sig = 0
        for contrast in self.contrasts:
            de_results = self.load_de_results(contrast)
            if de_results is not None:
                # Determine padj column name based on DE method
                padj_col = 'padj' if 'padj' in de_results.columns else 'FDR'
                fc_col = 'log2FoldChange' if 'log2FoldChange' in de_results.columns else 'logFC'
                
                # Count significant genes
                sig_genes = ((de_results[padj_col] < 0.05) & (abs(de_results[fc_col]) > 1)).sum()
                total_sig += sig_genes
        
        return total_sig
    
    def _get_sample_statistics(self):
        """Get alignment statistics for all samples."""
        sample_stats = {}
        
        for sample in self.samples:
            # Get flagstat info
            flagstat_file = os.path.join(self.result_path, "03_alignment", sample, f"{sample}_flagstat.txt")
            if os.path.exists(flagstat_file):
                with open(flagstat_file, 'r') as f:
                    stats = f.readlines()
                
                # Parse flagstat output
                total_reads = int(stats[0].split()[0])
                mapped_reads = int(stats[4].split()[0])
                mapping_rate = (mapped_reads / total_reads) * 100 if total_reads > 0 else 0
                
                sample_stats[sample] = {
                    'Total Reads': total_reads,
                    'Mapped Reads': mapped_reads,
                    'Mapping Rate (%)': round(mapping_rate, 2)
                }
        
        if sample_stats:
            return pd.DataFrame.from_dict(sample_stats, orient='index')
        
        return None
    
    def _load_sample_metadata(self):
        """Load sample metadata from config or metadata file."""
        # In a real implementation, this would parse the config file
        # or load a metadata file
        
        # For this example, we'll create mock metadata
        metadata = {}
        
        for sample in self.samples:
            # Extract condition from sample name (assuming naming convention)
            if "treatment" in sample.lower() or "treat" in sample.lower():
                condition = "treatment"
            elif "control" in sample.lower() or "ctrl" in sample.lower():
                condition = "control"
            else:
                # Default to sample name if no condition can be inferred
                condition = "sample"
            
            # Extract batch information (assuming naming convention)
            if "batch1" in sample.lower() or "b1" in sample.lower():
                batch = 1
            elif "batch2" in sample.lower() or "b2" in sample.lower():
                batch = 2
            else:
                batch = 1
            
            metadata[sample] = {
                'condition': condition,
                'batch': batch
            }
        
        if metadata:
            return pd.DataFrame.from_dict(metadata, orient='index')
        
        return None


# Main function to run the dashboard
def main():
    dashboard = RNASeqDashboard()
    dashboard.render_dashboard()


if __name__ == "__main__":
    main()
