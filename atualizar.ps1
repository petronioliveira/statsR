# Script PowerShell para atualizar Quarto Book no GitHub

Write-Host "🔄 Renderizando o livro..."
quarto render

Write-Host "📂 Adicionando arquivos modificados..."
git add _quarto.yml *.qmd referencias.bib styles/ docs

# Pergunta a mensagem do commit
$mensagem = Read-Host "📝 Digite a mensagem do commit"

Write-Host "📝 Criando commit com mensagem: $mensagem"
git commit -m "$mensagem"

Write-Host "🚀 Enviando para GitHub..."
git push origin main

Write-Host "✅ Atualização concluída! Confira em https://petronioliveira.github.io/statsR"